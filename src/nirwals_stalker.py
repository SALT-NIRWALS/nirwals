#!/usr/bin/env python3


import logging
import os
import argparse
import astropy.io.fits as pyfits

import multiparlog as mplog

import watchdog
from watchdog.observers import Observer
from watchdog.events import LoggingEventHandler, FileSystemEventHandler

import time
import nirwals

from nirwals import NIRWALS, dump_options



class NirwalsQuicklook(watchdog.events.PatternMatchingEventHandler):

    def __init__(self, holding_dir):
        self.logger = logging.getLogger("NirwalsQL")

        self.current_base = None
        self.zero_read = None
        self.good_sequence = False

        self.holding_dir = holding_dir

        # call the parent initialization routine
        # Set the patterns for PatternMatchingEventHandler
        watchdog.events.PatternMatchingEventHandler.__init__(
            self, patterns=['*.fits'],
            ignore_directories=True, case_sensitive=False)
        # super(NirwalsQuicklook, self).__init__()

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
        refpixels = nirwals.reference_pixels_to_background_correction(data)
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
        refpixels = nirwals.reference_pixels_to_background_correction(data)
        data -= refpixels

        # also subtract the minimum read
        data -= self.read_minimum

        # divide by exposure time
        exptime = hdulist[0].header['ACTEXP'] / 1000.
        data /= (exptime - self.read_mimimum_exptime)

        hdulist[0].data = data
        out_fn,ext = os.path.splitext(fn)
        out_fn += "__qred.fits"
        stage_fn = os.path.join(self.holding_dir, out_fn)
        self.logger.info("Writing quick-reduced frame to %s" % (stage_fn))
        hdulist.writeto(stage_fn, overwrite=True)

        pass

#        def on_created(self, event):
    def on_closed(self, event):
        print("Watchdog received closed event - % s." % event.src_path)
        # Event is created, you can process it now

        dir, fn = os.path.split(event.src_path)
        seq_base = ".".join(fn.split(".")[:-3])
        print(seq_base)
        print(type(event.src_path))
        if (seq_base == self.current_base):
            self.next_read(event.src_path)
        else:
            self.start_new_sequence(event.src_path)

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





if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsStalker")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--stage", dest="staging_dir", default="./",
                         help="staging directory")
    cmdline.add_argument("directory", nargs=1, help="name of directory to watch")
    args = cmdline.parse_args()

    path2watch = args.directory[0]
    print(path2watch)

    logging.info(f'start watching directory {path2watch!r}')
    event_handler = NirwalsQuicklook(args.staging_dir)
    observer = Observer()
    observer.schedule(event_handler, path2watch, recursive=False)
    observer.start()

    # time.sleep(1)
    # observer.start_new_sequence('test')
    # time.sleep(1)
    # observer.new_read("test2")

    try:
        while True:
            time.sleep(1)
    finally:
        observer.stop()
        observer.join()

