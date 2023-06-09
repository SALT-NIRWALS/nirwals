#!/usr/bin/env python3


import os
import sys
import numpy
import logging
import argparse

import astropy.io.fits as pyfits

import multiparlog as mplog

if __name__ == "__main__":

    mplog.setup_logging(debug_filename="../../debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("DarkCorrect")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--dark", dest="dark_fn", type=str, default=None,
                         help="calibration dark")
    cmdline.add_argument("--output", dest="output_postfix", type=str, default="reduced",
                         help="addition to output filename")

    cmdline.add_argument("--scaling", dest="scaling", default=None, type=float,
                         help="scaling factor")
    cmdline.add_argument("--dumpscaling", dest="dump_scaling",
                         default=False, action='store_true',
                         help="scaling factor")
    cmdline.add_argument("--fillbad", dest="fill_bad",
                         default=None, type=float,
                         help="fill bad pixels with fixed value")

    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()


    # in_fn = sys.argv[1]
    # dark_fn = sys.argv[2]
    # out_fn = sys.argv[3]

    logger.info("Using DARK %s" % (args.dark_fn))
    darkhdu = pyfits.open(args.dark_fn)
    dark = darkhdu['DARK.MASKED'].data

    for in_fn in args.files:

        logger.info("Working on %s" % (in_fn))
        hdulist = pyfits.open(in_fn)
        try:
            pers_signal = hdulist['PERS.SIGNAL'].data
        except KeyError:
            logger.error("Unable to find PERS.SIGNAL extension in %s, skipping" % (in_fn))
            continue

        corrected = pers_signal - dark
        if (args.fill_bad is not None):
            fill = ~numpy.isfinite(corrected) | (pers_signal < -9)
            corrected[fill] = args.fill_bad

        dir,bn = os.path.split(in_fn)
        bn_ext,ext = os.path.splitext(bn)
        out_fn = "%s_%s.fits" % (bn_ext, args.output_postfix)

        out_hdulist = pyfits.HDUList([
            hdulist[0],
            pyfits.ImageHDU(data=corrected, name='SCI.DARKCORR'),
        ])
        logger.info("Writing output to %s" % (out_fn))
        out_hdulist.writeto(out_fn, overwrite=True)

        if (args.dump_scaling):
            # assume all masked pixels
            masked_pixels = darkhdu['DARK.S2N'].data > 5
            pers = pers_signal[masked_pixels]
            _dark = dark[masked_pixels]
            # print(pers.shape, _dark.shape)

            dmp_fn = out_fn[:-5]+"__dark_factor.dmp"
            logger.info("writing scaling dump to %s" % (dmp_fn))
            numpy.savetxt(dmp_fn,
                          numpy.array([
                              _dark,
                              pers,
                              corrected[masked_pixels],
                              darkhdu['DARK.RMS'].data[masked_pixels],
                              darkhdu['DARK.S2N'].data[masked_pixels],
                          ]).T)

