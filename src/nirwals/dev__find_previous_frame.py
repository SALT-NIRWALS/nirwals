#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import logging
import multiparlog as mplog

# from rss_reduce import find_previous_frame
import nirwals_filepicker

if __name__ == "__main__":

    mplog.setup_logging(debug_filename="../../debug.log",
                        log_filename="run_analysis.log")

    sci_fn = sys.argv[1]
    dir_name = sys.argv[2]

    sci_hdu = pyfits.open(sci_fn)

    filepicker = rss_filepicker.PreviousFilePicker(dir_name)
    # result = find_previous_frame(sci_hdu[0].header, dir_name)
    fn, delta_t = filepicker.find_closest(sci_hdu[0].header)

    print(fn)

