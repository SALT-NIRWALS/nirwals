#!/usr/bin/env python3


import logging
import multiprocessing
import os
import sys
import argparse
import astropy.io.fits as pyfits
import astropy.samp.errors
import pyvo.samp as sampy

import multiparlog as mplog
import numpy

import watchdog
import watchdog.events
import watchdog.observers

import sysv_ipc
import time
import nirwals


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

    def on_closed(self, event):
        print("Watchdog received closed event - % s." % event.src_path)
        # Event is created, you can process it now

        # Hand on work to parallel worker
        self.handler_queue.put(event.src_path)


class NirwalsOnTheFlyReduction(multiprocessing.Process):

    def __init__(self, incoming_queue, staging_dir, nonlinearity_file, refpixelmode, samp_cli=None, shmem_ds9=None, force_write=False):
        super(NirwalsOnTheFlyReduction, self).__init__()

        self.logger = logging.getLogger("WatchdogProcess")

        self.current_base = None
        self.zero_read = None
        self.good_sequence = False
        self.sequence_firstread_fn = None
        self.read_minimum = None
        self.read_minimum_exptime = None
        self.latest_result = None
        self.saturated = None
        self.force_write = force_write

        self.staging_dir = staging_dir
        self.incoming_queue = incoming_queue
        self.refpixelmode = refpixelmode
        self.latest_result = None

        self.nonlin_poly = None
        self.nonlin_poly_order = 0
        self.nonlinearity_fn = nonlinearity_file
        if (nonlinearity_file is not None and os.path.isfile(nonlinearity_file)):
            hdulist = pyfits.open(nonlinearity_file)
            self.nonlin_poly = hdulist[0].data
            self.nonlin_poly_order = self.nonlin_poly.shape[0] - 1
            self.logger.debug("Read non-linearity corrections from %s (order: %d)" % (
                nonlinearity_file, self.nonlin_poly_order
            ))

        self.samp_cli = samp_cli

        self.shmem_ds9 = shmem_ds9
        self.shmem_buf = None
        if (shmem_ds9 is not None):
            class Dummy(object): pass
            d = Dummy()
            d.__array_interface__ = {
                 'data' : (shmem_ds9.address, False),
                 'typestr' : "=f4", #FloatType, #"uint8", #numpy.uint8.str,
                 'descr' : "", #"UINT8", #numpy.uint8.descr,
                 'shape' : (2048,2048), #(n_bytes//4,), #(self.n_bytes/4), #
                 'strides' : None,
                 'version' : 3
            }
            self.shmem_buf = numpy.asarray(d) #.reshape((ny,nx))
            if (samp_cli is not None):
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='frame 7')
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='shm array shmid %d [xdim=%d,ydim=%d,bitpix=-32]' % (
                    self.shmem_ds9.id,2048,2048))

    def start_new_sequence(self, filename):
        self.logger.info("Starting new sequence: %s" % (filename))

        dirname, fn = os.path.split(filename)
        seq_base = ".".join(fn.split(".")[:-3])

        # update the name of the sequence we are currently working on
        self.current_base = seq_base

        # find the first read in this sequence
        self.sequence_firstread_fn = os.path.join(dirname, "%s.1.1.fits" % (seq_base))
        if (not os.path.isfile(self.sequence_firstread_fn)):
            self.logger.warning("Unable to find first read in this sequence (%s)" % (self.sequence_firstread_fn))
            self.good_sequence = False
            return

        hdulist = pyfits.open(self.sequence_firstread_fn)
        data = hdulist[0].data.astype(float)
        refpixels = nirwals.reference_pixels_to_background_correction(data, mode=self.refpixelmode)
        data -= refpixels

        self.read_minimum = data
        self.read_mimimum_exptime = hdulist[0].header['ACTEXP'] / 1e6
        self.logger.info("Sequence start is GOOD!")
        self.good_sequence = True

        self.latest_result = None
        self.saturated = None
        pass

    def next_read(self, filename):
        self.logger.debug("Handling new read: %s (%s)" % (filename, type(filename)))

        _, fn = os.path.split(filename)

        hdulist = pyfits.open(filename)
        if ('WATCHRED' in hdulist[0].header):
            if (hdulist[0].header['WATCHRED'] == 'yes'):
                self.logger.warning("Found watchdog output file among the input files. Check configuration!!!!")
                return None

        data = hdulist[0].data.astype(float)
        saturated = (data > 62000)
        if (self.saturated is None):
            self.saturated = saturated
        else:
            self.saturated = saturated | self.saturated
            # make sure that pixels that were marked at saturated in any read stay marked as saturated in all
            # subsequent reads as well

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
            self.logger.debug("Done with non-linearity correction (%.4f)" % (t2-t1))

        # divide by exposure time
        exptime = hdulist[0].header['ACTEXP'] / 1e6
        data /= (exptime - self.read_mimimum_exptime)

        if (self.latest_result is None):
            self.latest_result = data
        else:
            self.latest_result[~self.saturated] = data[~self.saturated]

        hdulist[0].data = self.latest_result

        # Put our special mark on this frame to ensure we don't accidentally process the output of this file again
        hdulist[0].header['WATCHRED'] = ("yes", "file was processed using the nirwals watchdog")

        out_fn = None
        if (self.shmem_ds9 is None or self.force_write):
            out_fn, ext = os.path.splitext(fn)
            out_fn += "__qred.fits"
            stage_fn = os.path.abspath(os.path.join(self.staging_dir, out_fn))
            self.logger.debug("Writing quick-reduced frame to %s" % (stage_fn))
            hdulist.writeto(stage_fn, overwrite=True)

        if (self.samp_cli is not None):
            # samp_msg = {
            #     'samp.mtype': 'image.load.fits',
            #     'samp.params': {
            #         'url': 'file://'+os.path.abspath(stage_fn),
            #         'name': 'latest NIRWALS frame',
            #     }
            # }
            # self.samp_cli.notify_all(samp_msg)
            if (self.shmem_ds9 is not None):
                self.shmem_buf[:,:] = data[:,:]
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='update now')
                self.logger.debug("Using shared memory, telling ds9 to update")
            else:
                self.logger.debug("Sending command to load new FITS (%s) to ds9" % (stage_fn))
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='frame 7')
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='preserve pan yes')
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='preserve region yes')
                self.samp_cli.enotify_all(mtype='ds9.set', cmd='fits %s' % (stage_fn))
                # self.samp_cli.enotify_all(mtype='ds9.set', cmd='array file://%s' % (os.path.abspath(stage_fn)))

        return out_fn

    def run(self):

        while (True):
            t1 = time.time()
            try:
                job = self.incoming_queue.get()
            except KeyboardInterrupt:
                self.logger.info("Shutting down on-the-fly worker after user command")
                break
            except Exception as e:
                self.logger.critical(str(e))
                self.incoming_queue.task_done()
                break
            t2 = time.time()
            self.logger.debug("got new task after waiting for %.3f seconds" % (t2-t1))

            if (job is None):
                self.logger.info("Shutting down worker")
                self.incoming_queue.task_done()
                break

            start_time = time.time()
            new_filename = job
            _, fn = os.path.split(new_filename)
            seq_base = ".".join(fn.split(".")[:-3])
            # print(seq_base)
            # print(type(new_filename))
            if (seq_base == self.current_base):
                out_fn = self.next_read(new_filename)
            else:
                self.start_new_sequence(new_filename)
                out_fn = "NEW SEQ"
            end_time = time.time()

            self.incoming_queue.task_done()
            self.logger.info("On-the-fly processing of %s --> %s completed after %.3f seconds" % (
                new_filename, out_fn, end_time-start_time))


samp_metadata = {
    "samp.name": "NIRWALS_watchdog",
    "samp.description.text": "NIRWALS directory watcher and on-the-fly reduction",
    "samp.icon.url": "https://avatars.githubusercontent.com/u/95205451?s=48&v=4",
    "samp.documentation.url": "https://github.com/SALT-NIRWALS/nirwals/",
    "author.name": "Ralf Kotulla",
    "author.email": "kotulla@wisc.edu",
    "author.affiliation": "University of Wisconsin - Madison",
    "home.page": "https://github.com/rkotulla",
    "cli1.version": "0.01",
}


def main():

    mplog.setup_logging(debug_filename="nirwals_debug.log",
                        log_filename="nirwals_watchdog.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsWatchdog")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--stage", dest="staging_dir", default="./",
                         help="staging directory")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--refpixel", dest="ref_pixel_mode", default='blockyslope2',
                         help="reference pixels mode [default: NO]")
    cmdline.add_argument("--test", dest="test", default=None,
                         help="test mode; syntax: --test=delay:@filelist")
    cmdline.add_argument("--nowait", dest="no_wait_for_samp", default=False, action='store_true',
                         help="do not wait for SAMP server")
    cmdline.add_argument("--shmem", dest="use_shared_memory", default=False, action='store_true',
                         help="use shared memory to display files in ds9")
    cmdline.add_argument("--write", dest="write_reduced", default=False, action='store_true',
                         help="write reduced frames (even if not needed)")
    cmdline.add_argument("directory", nargs=1, help="name of directory to watch")
    args = cmdline.parse_args()

    path2watch = os.path.abspath(args.directory[0])

    stage_dir = os.path.abspath(args.staging_dir)
    if (not os.path.isdir(path2watch)):
        logger.error("Directory to watch (%s) not found, please check & rerun" % (path2watch))
        sys.exit(0)
    if (not os.path.isdir(stage_dir)):
        logger.error("Staging directory (%s) does not exist, please create it first and re-run" % (stage_dir))
        sys.exit(0)
    if (os.path.samefile(path2watch, stage_dir)):
        logger.error("Watch-directory and staging directory are identical, this should not happen")
        logger.error("Please adjust watch directory and/or staging directory [via the --stage option]")
        sys.exit(0)

    # print(path2watch)

    logger.info("Connecting to ds9 via SAMP protocol")
    samp_cli = None
    start_waiting = time.time()
    last_update = 0
    while (True):
        try:
            samp_cli = sampy.SAMPIntegratedClient(metadata=samp_metadata)
            samp_cli.connect()
            if (not samp_cli.is_connected):
                samp_cli = None
        except astropy.samp.errors.SAMPHubError as e:
            samp_cli = None
        except Exception as e:
            logger.critical("Error while establishing SAMP link: %s" % (str(e)))
            pass
        if (samp_cli is None):
            if (args.no_wait_for_samp):
                logger.critical("No SAMPHub found, automatic forwarding to ds9 disabled")
                break
            else:
                t = time.time()
                if (last_update == 0):
                    logger.info("No SAMPHub found, waiting for SAMPhub to come online")
                    last_update = t
                elif ((t-last_update) > 60):
                    logger.info("Still no SAMPHub found, continuing to wait for SAMPhub to come online (total %d seconds by now)" % (t-start_waiting))
                    last_update = t
                time.sleep(1)
        else:
            logger.info("Successfully connected to SAMPhub")
            break

    # Open a frame in ds9
    if (samp_cli is not None):
        samp_cli.enotify_all(mtype='ds9.set', cmd='frame 7')

    shmem_ds9 = None
    if (args.use_shared_memory):
        logger.info("Allocating shared memory for interaction with ds9")
        dummy = numpy.array([], dtype=numpy.float32)
        nx, ny = 2048, 2048
        n_bytes = dummy.itemsize * nx * ny
        shmem_ds9 = sysv_ipc.SharedMemory(
            key=None, mode=0o666, size=n_bytes,
            flags=sysv_ipc.IPC_CREAT | sysv_ipc.IPC_EXCL)
    else:
        logger.info("Writing reduced files to forward to ds9")

    logger.info("Starting on-the-fly reduction process")
    job_queue = multiprocessing.JoinableQueue()
    watchdog_worker = NirwalsOnTheFlyReduction(
        job_queue, staging_dir=stage_dir,
        nonlinearity_file=args.nonlinearity_fn,
        refpixelmode=args.ref_pixel_mode,
        samp_cli=samp_cli,
        shmem_ds9=shmem_ds9,
        force_write=args.write_reduced,
    )
    watchdog_worker.daemon = True
    watchdog_worker.start()

    logger.info("Starting to watch directory for new files: %s" % (path2watch))
    event_handler = NirwalsQuicklook(job_queue)
    observer = watchdog.observers.Observer()
    observer.schedule(event_handler, path2watch, recursive=False)
    observer.start()

    # For testing, simulate the arrival of a bunch of new files
    if (args.test is not None):
        testdelay = float(args.test.split(":")[0])
        list_fn = args.test.split("@")[1]
        testfiles = []
        with open(list_fn, "r") as lf:
            lines = lf.readlines()
            for l in lines:
                if (l.strip().startswith("#") or len(l.strip()) < 1):
                    continue
                fn = l.strip().split()[0]
                if (os.path.isfile(fn)):
                    testfiles.append(fn)
        logger.info("Test-mode activated, feeding %d files with delays of %.2f seconds" % (len(testfiles), testdelay))
        for fn in testfiles:
            time.sleep(testdelay)
            job_queue.put(fn)

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        logger.debug("Starting shutdown proceduce")
        pass
    finally:
        # Shut down the directory watchdog
        logger.debug("Stopping watchdog")
        observer.stop()
        observer.join()

        # Terminate the on-the-fly reduction worker
        logger.debug("Shutting down on-the-fly reduction")
        job_queue.put(None)
        job_queue.close()
        watchdog_worker.join()

    logger.info("All shut down, good bye & have a great day!")

if __name__ == "__main__":
    main()