#!/usr/bin/env python3

import sys
import os
import multiparlog as mplog
import argparse
import logging

# for verification
import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    cmdline.add_argument("--refpixel", dest="ref_pixel_mode", default='blockyslope2',
                         help="reference pixels mode [default: NO]")
    cmdline.add_argument("--verify", dest="verify", default=False, action='store_true',
                         help="verify results rather than fitting coefficients")
    cmdline.add_argument("--pixels", dest="pixels", type=str,
                         help="list of pixel coordinates or file with coordinates")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    fn = args.files[0]
    saturation_fn = args.saturation
    rss = NIRWALS(fn, saturation=saturation_fn,
                  max_number_files=args.max_number_files,
                  use_reference_pixels=args.ref_pixel_mode,)

    # rss.reduce(write_dumps=False)
    # rss.write_results()

    rss.load_all_files()
    # rss.subtract_first_read()

    if (not args.verify):
        rss.fit_nonlinearity(ref_frame_id=4, make_plot=False)

    else:
        rss.reduce(dark_fn=None,)
        # Verify rather than fit results
        coeff_hdu = pyfits.open(args.nonlinearity_fn)
        coeffs = coeff_hdu[0].data
        print(coeffs.shape)

        # Read all pixel coordinates, either from command line or @file
        pixels = []
        if (args.pixels.startswith("@")):
            with open(args.pixels[1:]) as f:
                lines = f.readlines()
                for l in lines:
                    if (l.strip().startswith("#")):
                        continue
                    items = l.split()
                    xy = [int(round(float(x))) for x in items[0:2]]
                    pixels.append(xy)
        else:
            pairs = args.pixels.split(":")
            for p in pairs:
                items = p.split(",")
                xy = [int(round(float(x))) for x in items[0:2]]
                pixels.append(xy)
        print(pixels)

        with PdfPages("nonlinearity_verification.pdf") as pdf:
            for xy in pixels:
                [x,y] = xy
                print(xy,x,y)

                fig = plt.figure()
                ax = fig.add_subplot(111)

                raw_sequence = rss.linearized_cube[:,y,x]
                read_number = numpy.arange(raw_sequence.shape[0])
                print(raw_sequence.shape)
                ax.scatter(read_number, raw_sequence, s=3, alpha=0.2)
                ax.plot(read_number, raw_sequence, 'b-', linewidth=1)

                fig.suptitle("Pixel %d , %d" % (x,y))

                pdf.savefig(fig)
                # fig.savefig()
                #

    # rss.plot_pixel_curve(818,1033)
    # rss.plot_pixel_curve(1700,555)
