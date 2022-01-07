#!/usr/bin/env python3

import sys
import os
import argparse

import astropy.io.fits as pyfits

import rss_reduce


if __name__ == "__main__":

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--output", dest="output_fn", type=str, default=None,
                         help="output filename")
    cmdline.add_argument("--dumps", dest="write_dumps", default=False, action='store_true',
                         help="write intermediate process data [default: NO]")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    for fn in args.files:

        rss = rss_reduce.RSS(fn, max_number_files=args.max_number_files)
        if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
            rss.read_nonlinearity_corrections(args.nonlinearity_fn)
        rss.reduce(write_dumps=args.write_dumps)

        # Figure out what the incremental exposure time per read is
        exptime = rss.first_header['USEREXP'] / 1000. # raw values are in milli-seconds
        delta_exptime = exptime / rss.first_header['NGROUPS']

        darkrate = rss.weighted_mean / delta_exptime

        dark_hdu = pyfits.HDUList([
            pyfits.PrimaryHDU(header=rss.first_header),
            pyfits.ImageHDU(data=darkrate, name='DARKRATE')
        ])

        if args.output_fn is None:
            out_fn = rss.filebase + ".darkrate.fits"
        else:
            out_fn = args.output_fn
        print("Writing darkrate image to %s ..." % (out_fn))
        dark_hdu.writeto(out_fn, overwrite=True)

