#!/usr/bin/env python3

import sys
import os
import argparse
import numpy
import scipy
import scipy.optimize
import matplotlib.pyplot as plt

import astropy.io.fits as pyfits

import nirwals_reduce

if __name__ == "__main__":

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--output", dest="output_fn", type=str, default=None,
                         help="output filename")
    cmdline.add_argument("--op", dest="operation", type=str, default='-',
                         help="output filename")

    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    if (len(args.files) != 2):
        print("Need exactly two input files, quitting")
        sys.exit(0)

    datacubes = []
    for fn in args.files:
        print("Reading input %s" % (fn))
        rss = rss_reduce.RSS(fn, max_number_files=args.max_number_files)
        rss.load_all_files()
        rss.reset_frame = rss.image_stack[0]
        if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
            rss.read_nonlinearity_corrections(args.nonlinearity_fn)

        reset_frame_subtracted = rss.image_stack - rss.reset_frame
        # apply any necessary corrections for nonlinearity and other things
        print("Applying non-linearity corrections")
        linearized = rss.apply_nonlinearity_corrections(reset_frame_subtracted)
        # print("linearized = ", linearized)
        if (linearized is None):
            print("No linearized data found, using raw data instead")
            linearized = reset_frame_subtracted

        datacubes.append(linearized)

    outcube = None
    if (args.operation == "-" or args.operation.lower()=='minus'):
        outcube = datacubes[0] - datacubes[1]
    elif (args.operation == "+"):
        outcube = datacubes[0] + datacubes[1]
    elif (args.operation == "x"):
        outcube = datacubes[0] * datacubes[1]
    elif (args.operation == "/"):
        outcube = datacubes[0] / datacubes[1]
    else:
        print("operation (%s) not understood" % (args.operation))

    if (outcube is not None):
        out_fn = "cubemath_output.fits" if args.output_fn is None else args.output_fn
        print("Writing results to %s" % (out_fn))
        pyfits.PrimaryHDU(data=outcube).writeto(out_fn, overwrite=True)

    print("all done!")