#!/usr/bin/env python3


import logging
import multiprocessing
import os
import argparse
import astropy.io.fits as pyfits

import multiparlog as mplog
import numpy

import watchdog
from watchdog.observers import Observer
from watchdog.events import LoggingEventHandler, FileSystemEventHandler

import time
import nirwals

from nirwals import NIRWALS, dump_options



class NirwalsQuicklook(watchdog.events.PatternMatchingEventHandler):

    def __init__(self, handler_queue):
        self.logger = logging.getLogger("NirwalsQL")
        self.handler_queue = handler_queue

        # call the parent initialization routine
        # Set the patterns for PatternMatchingEventHandler
        watchdog.events.PatternMatchingEventHandler.__init__(
            self, patterns=['*.fits'],
            ignore_directories=True, case_sensitive=False)
        super(NirwalsQuicklook, self).__init__()


#        def on_created(self, event):
    def on_closed(self, event):
        print("Watchdog received closed event - % s." % event.src_path)
        # Event is created, you can process it now

        # Hand on work to parallel worker
        self.handler_queue.put(event.src_path)


        # return event
            # this is a new sequence

    # @staticmethod
    # def on_any_event(event):
    #     if event.is_directory:
    #         return None
    #     elif event.event_type == 'created':
    #         # Event is created, you can process it now
    #         print("Watchdog received created event - % s." % event.src_path)
    #     else:
    #         return None
    #
    #     if (not event.src_path.endswith(".fits")):
    #         return None



class NirwalsOnTheFlyReduction(multiprocessing.Process):

    def __init__(self, incoming_queue, staging_dir, nonlinearity_file, refpixelmode, ds9):
        super(NirwalsOnTheFlyReduction, self).__init__()

        self.logger = logging.getLogger("StalkerProcess")

        self.current_base = None
        self.zero_read = None
        self.good_sequence = False

        self.staging_dir = staging_dir
        self.incoming_queue = incoming_queue
        self.refpixelmode = refpixelmode

        self.nonlin_poly = None
        self.nonlin_poly_order = 0
        self.nonlinearity_fn = nonlinearity_file
        if (nonlinearity_file is not None and os.path.isfile(nonlinearity_file)):
            hdulist = pyfits.open(nonlinearity_file)
            self.nonlin_poly = hdulist[0].data
            self.nonlin_poly_order = self.nonlin_poly.shape[0] - 1
            self.logger.info("Read non-linearity corrections from %s (order: %d)" % (
                nonlinearity_file, self.nonlin_poly_order
            ))

    def start_new_sequence(self, filename):
        print("Starting new sequence: %s" % (filename))

        dir, fn = os.path.split(filename)
        seq_base = ".".join(fn.split(".")[:-3])

        # update the new of the sequence we are currently working on
        self.current_base = seq_base

        # find the first read in this sequence
        self.sequence_firstread_fn = os.path.join(dir, "%s.1.1.fits" % (seq_base))
        if (not os.path.isfile(self.sequence_firstread_fn)):
            self.logger.warning("Unable to find first read in this sequence (%s)" % (self.sequence_firstread_fn))
            self.good_sequence = False
            return

        hdulist = pyfits.open(self.sequence_firstread_fn)
        data = hdulist[0].data.astype(float)
        refpixels = nirwals.reference_pixels_to_background_correction(data, mode=self.refpixelmode)
        data -= refpixels

        self.read_minimum = data
        self.read_mimimum_exptime = hdulist[0].header['ACTEXP'] / 1000.
        self.logger.info("Sequence start is GOOD!")
        self.good_sequence = True
        pass

    def next_read(self, filename):
        print("Handling new read: %s (%s)" % (filename, type(filename)))

        dir, fn = os.path.split(filename)

        hdulist = pyfits.open(filename)
        data = hdulist[0].data.astype(float)
        refpixels = nirwals.reference_pixels_to_background_correction(data, mode=self.refpixelmode)
        data -= refpixels

        # also subtract the minimum read
        data -= self.read_minimum

        # apply nonlinearity correction
        if (self.nonlin_poly is not None):
            t1 = time.time()
            nlc = numpy.zeros_like(data, dtype=float)
            for p in range(self.nonlin_poly_order):
                nlc += self.nonlin_poly[p] * numpy.power(data, self.nonlin_poly_order-p)
            data = nlc
            t2 = time.time()
            self.logger.info("Done with non-linearity correction (%.4f)" % (t2-t1))

        # divide by exposure time
        exptime = hdulist[0].header['ACTEXP'] / 1000.
        data /= (exptime - self.read_mimimum_exptime)

        hdulist[0].data = data
        out_fn,ext = os.path.splitext(fn)
        out_fn += "__qred.fits"
        stage_fn = os.path.join(self.staging_dir, out_fn)
        self.logger.info("Writing quick-reduced frame to %s" % (stage_fn))
        hdulist.writeto(stage_fn, overwrite=True)

        pass

    def run(self):

        while (True):
            try:
                job = self.incoming_queue.get()
            except Exception as e:
                self.logger.critical(str(e))

            if (job is None):
                self.logger.info("Shutting down worker")
                self.incoming_queue.task_done()
                break

            start_time = time.time()
            new_filename = job
            dir, fn = os.path.split(new_filename)
            seq_base = ".".join(fn.split(".")[:-3])
            print(seq_base)
            print(type(new_filename))
            if (seq_base == self.current_base):
                self.next_read(new_filename)
            else:
                self.start_new_sequence(new_filename)
            end_time = time.time()

            self.incoming_queue.task_done()
            self.logger.info("On-the-fly processing of %s completed after %.3f seconds" % (
                new_filename, end_time-start_time))

if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsStalker")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--stage", dest="staging_dir", default="./",
                         help="staging directory")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--refpixel", dest="ref_pixel_mode", default='blockyslope2',
                         help="reference pixels mode [default: NO]")
    cmdline.add_argument("directory", nargs=1, help="name of directory to watch")
    args = cmdline.parse_args()

    path2watch = args.directory[0]
    print(path2watch)

    job_queue = multiprocessing.JoinableQueue()
    stalker_worker = NirwalsOnTheFlyReduction(
        job_queue, staging_dir=args.staging_dir,
        nonlinearity_file=args.nonlinearity_fn,
        refpixelmode=args.ref_pixel_mode,
        ds9=None)
    stalker_worker.daemon = True
    stalker_worker.start()

    logging.info(f'start watching directory {path2watch!r}')
    event_handler = NirwalsQuicklook(job_queue)
    observer = Observer()
    observer.schedule(event_handler, path2watch, recursive=False)
    observer.start()

    for n in range(6):
        time.sleep(00.7)
        job_queue.put('/nas/t7black/salt/incoming/N202303080003.3.1.%d.fits' % (n+1))
    # time.sleep(1)
    # job_queue.put('/nas/t7black/salt/incoming/N202303080003.3.1.2.fits')
    # time.sleep(1)
    # job_queue.put('/nas/t7black/salt/incoming/N202303080003.3.1.3.fits')
    # observer.new_read("test2")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        pass
    finally:
        # Shut down the directory watchdog
        print("Stopping watchdog")
        observer.stop()
        observer.join()

        # Terminate the on-the-fly reduction worker
        print("Shutting down on-the-fly reduction")
        job_queue.put(None)
        job_queue.close()
        stalker_worker.join()

