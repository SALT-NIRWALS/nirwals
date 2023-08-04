#!/usr/bin/env python3

import sys
# print(sys.path)

import logging
import os

import multiparlog as mplog

import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import itertools
import multiprocessing
import multiprocessing.shared_memory
import argparse

from astropy import log
log.setLevel('ERROR')

import warnings
warnings.filterwarnings('ignore')

import astropy
print(astropy.__path__)

from nirwals import NIRWALS, dump_options

if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsReduce")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--flat", dest="flatfield_fn", type=str, default=None,
                         help="calibration flatfield")
    cmdline.add_argument("--dark", dest="dark_fn", type=str, default=None,
                         help="calibration dark")
    cmdline.add_argument("--output", dest="output_postfix", type=str, default="reduced",
                         help="addition to output filename")

    cmdline.add_argument("--persistency", dest="persistency_mode", type=str, default="quick",
                         help="persistency mode")
    cmdline.add_argument("--saturation", dest="saturation", default=62000,
                         help="saturation value/file")

    # cmdline.add_argument("healpix32", nargs="+", type=int,
    #                      help="list of input filenames")
    # cmdline.add_argument("--rerun", type=int, default=6,
    #                      help="rerun")
    cmdline.add_argument("--dumps", dest="write_dumps", default=None,
                         help="write intermediate process data [default: NO]")
    cmdline.add_argument("--debugpngs", dest="write_debug_pngs", default=False, action='store_true',
                         help="generate debug plots for all pixels with persistency [default: NO]")
    cmdline.add_argument("--refpixel", dest="ref_pixel_mode", default='none',
                         help="reference pixels mode [default: NO]")
    cmdline.add_argument("--flat4salt", dest="write_flat_for_salt", default=False, action='store_true',
                         help="write a flat, 1-extension FITS file for SALT")
    cmdline.add_argument("--report", dest="report_provenance", default=False, action='store_true',
                         help="report ata provenance at end of processing")
    cmdline.add_argument("--speedy", dest="speedy", default=False, action='store_true',
                         help="speed up processing by adaptively reducing sample reads")
    cmdline.add_argument("--ncores", dest="n_cores", default=None, type=int,
                         help="number of CPU cores to use")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    dumpfiles = []
    if (args.write_dumps is not None):
        for di in [x.lower() for x in args.write_dumps.split(",")]:
            if (di in dump_options):
                dumpfiles.append(di)

    for fn in args.files:
        # fn = sys.argv[1]

        try:
            rss = NIRWALS(fn, max_number_files=args.max_number_files,
                      use_reference_pixels=args.ref_pixel_mode,
                      saturation=args.saturation,
                      nonlinearity=args.nonlinearity_fn,
                      logger_name="Nirwals",
                      speedy=args.speedy,
                      n_cores=args.n_cores,
                      dumps=dumpfiles,
                      )
        except:
            logger.critical("Unable to start processing, read and resolve error message before continuing")
            continue

        # if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
        #     logger.info("Attempting to load non-linearity from %s" % (args.nonlinearity_fn))
        #     rss.read_nonlinearity_corrections(args.nonlinearity_fn)
        rss.reduce(dark_fn=args.dark_fn,
                   )

        # persistency_options = args.persistency_mode.split(":")
        # persistency_mode = persistency_options[0].lower()
        # have_persistency_results = False
        # if (persistency_mode == "none"):
        #     logger.info("Nothing to do")
        # elif (persistency_mode == "full"):
        #     logger.info("Calculating persistency results for all pixels")
        #     rss.fit_signal_with_persistency(previous_frame=None)
        #     have_persistency_results = True
        # elif (persistency_mode == "best"):
        #     logger.info("Using on-demand persistency fitting")
        #     if(len(persistency_options) < 2):
        #         logger.info("Insufficient information, defaulting to running on all pixels")
        #         rss.fit_signal_with_persistency(previous_frame=None)
        #         have_persistency_results = True
        #     else:
        #         opt = persistency_options[1]
        #         if (os.path.isfile(opt)):
        #             logger.info("Using optimized persistency mode (ref-fn: %s)" % (opt))
        #             rss.fit_signal_with_persistency(previous_frame=opt)
        #             have_persistency_results = True
        #         elif (os.path.isdir(opt)):
        #             logger.info("Searching for optimal reference frame in --> %s <--" % (opt))
        #             xxx = rss.find_previous_exposure(opt)  #find_previous_frame(rss.ref_header, opt)
        #             # print(xxx)
        #             ref_fn, delta_t = xxx #rss.find_previous_exposure(opt)  #find_previous_frame(rss.ref_header, opt)
        #             if (ref_fn is not None):
        #                 logger.info("Using optimized persistency mode using automatic ref-fn: %s (Dt=%.3f)" % (ref_fn, delta_t))
        #                 rss.fit_signal_with_persistency(previous_frame=ref_fn,
        #                                                 write_test_plots=args.write_debug_pngs)
        #                 have_persistency_results = True
        #             else:
        #                 logger.warning("No previous frame found, skipping persistency fit")
        #                 # rss.fit_signal_with_persistency(previous_frame=ref_fn)
        #                 # have_persistency_results = True
        #         else:
        #             logger.info("Unknown option to best mode (found: %s), skipping persistency modeling" % (opt))
        # else:
        #     logger.info("Unknown persistency request (%s)" % (persistency_mode))
        #
        # if (have_persistency_results):
        #     out_tmp = pyfits.PrimaryHDU(data=rss.persistency_fit_global)
        #     fit_fn = "%s.%s.persistencyfit.fits" % (rss.filebase, args.output_postfix)
        #     logger.info("Writing persistency fit to %s ..." % (fit_fn))
        #     out_tmp.writeto(fit_fn, overwrite=True)

        red_fn = "%s.%s.fits" % (rss.filebase, args.output_postfix)
        logger.info("Writing reduction results to %s" % (os.path.abspath(red_fn)))
        rss.write_results(fn=red_fn, flat4salt=args.write_flat_for_salt)

        if (args.report_provenance):
            rss.provenance.report()

        # rss.plot_pixel_curve(818,1033)
        # rss.plot_pixel_curve(1700,555)
        # rss.plot_pixel_curve(505,1660)

        del rss
        logger.info("all done!")
