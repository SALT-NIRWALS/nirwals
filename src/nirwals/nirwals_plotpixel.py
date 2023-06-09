#!/usr/bin/env python3

import sys
import os
import numpy
import astropy.io.fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
import argparse

import nirwals_reduce


if __name__ == "__main__":

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--output", dest="output_fn", type=str, default=None,
                         help="output filename")
    cmdline.add_argument("--plotdir", dest="plot_directory", type=str, default=None,
                         help="output filename")
    cmdline.add_argument("--dumps", dest="write_dumps", default=False, action='store_true',
                         help="write intermediate process data [default: NO]")
    cmdline.add_argument("--cumulative", dest="plot_cumulative", default=False, action='store_true',
                         help="make cumulative pixel plot [default: NO]")
    cmdline.add_argument("--diff", dest="plot_differential", default=False, action='store_true',
                         help="make differential pixel plot [default: NO]")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    cmdline.add_argument("--pixels", dest='pixels', nargs="+",
                         help="list of pixel coordinates")
    args = cmdline.parse_args()

    if (args.write_dumps):
        print("File-dumping enabled!")

    for fn in args.files:

        rss = rss_reduce.RSS(fn)
        if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
            rss.read_nonlinearity_corrections(args.nonlinearity_fn)
        rss.reduce(write_dumps=args.write_dumps,
                   dark_fn=None)
        # rss.write_results(fn="%s.%s.fits" % (rss.filebase, args.output_postfix))

        for pixel in args.pixels:
            print("Making plot for %s" % (pixel))
            xy = [int(s) for s in pixel.split(",")]
            if (len(xy) != 2):
                print("Can't understand pixel position: %s" % (pixel))
                continue
            else:
                if (args.plot_directory is not None):
                    filebase = os.path.join(args.plot_directory, rss.filebase)
                else:
                    filebase - rss.filebase
                rss.plot_pixel_curve(xy[0], xy[1], filebase=filebase,
                                     cumulative=args.plot_cumulative,
                                     differential=args.plot_differential,
                                     diff_vs_cum=True, show_errors=True)

