#!/usr/bin/env python3

import sys
import os
import multiparlog as mplog
import argparse
import logging

from nirwals import NIRWALS


if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsFitNonlinearity")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--saturation", dest="saturation", default=62000,
                         help="saturation value/file")
    cmdline.add_argument("--refpixel", dest="use_ref_pixels", default=False, action='store_true',
                         help="use reference pixels [default: NO]")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    fn = args.files[0]
    saturation_fn = args.saturation
    rss = NIRWALS(fn, saturation=saturation_fn,
                  max_number_files=args.max_number_files)

    # rss.reduce(write_dumps=False)
    # rss.write_results()

    rss.load_all_files()
    # rss.subtract_first_read()

    rss.fit_nonlinearity(ref_frame_id=4, make_plot=False)

    # rss.plot_pixel_curve(818,1033)
    # rss.plot_pixel_curve(1700,555)
