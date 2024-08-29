#!/usr/bin/env python3

import numpy
import scipy
import nirwals
import os
import sys
import matplotlib.pyplot as plt
import time
import astropy.io.fits as pyfits
import pandas
import multiparlog as mplog
import logging
import argparse

from astropy import log
log.setLevel('ERROR')

if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsReduce")

    files = ['linearized_datacube','nonlinearity_corrections','raw_datacube','results_datacube']
    for f in files:
        fn = os.path.join('/dev/shm/', f)
        if (os.path.isfile(fn)):
            print("deleting", fn)
            os.remove(fn)

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--output", dest="output_fn", type=str, default="eff_noise_prep.fits",
                         help="addition to output filename")
    cmdline.add_argument("--refpixel", dest="ref_pixel_mode", default='none',
                         help="reference pixels mode [default: NO]")
    cmdline.add_argument("--ncores", dest="n_cores", default=None, type=int,
                         help="number of CPU cores to use")
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--algorithm", dest="algorithm", default='linreg', type=str,
                         choices={"linreg", "rauscher2007", 'pairwise', 'linreg+recombine', 'rauscher2007+recombine'},
                         help="number of CPU cores to use")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()


    rss = nirwals.NIRWALS(
        fn=args.files[0],
        max_number_files=args.max_number_files,
        use_reference_pixels=args.ref_pixel_mode,
        saturation=53000,
        nonlinearity=args.nonlinearity_fn,
        logger_name="Nirwals",
        speedy=False,
        n_cores=args.n_cores,
        dumps=None,
    )

    # prepare all data
    rss.load_all_files()
    rss.apply_reference_pixel_corrections()
    rss.apply_nonlinearity_corrections()

    # start up output file
    out_hdulist = [pyfits.PrimaryHDU()]

    n = numpy.arange(1,20,1)
    n = numpy.append(n, numpy.arange(20,50,3))
    n = numpy.append(n, numpy.arange(50,824,5))
    n_groups = n

    logger.info("Preparing analysis, output will be in %s" % (args.output_fn))

    for grps in n_groups:
        print("\nWorking on %d groups\n" % (grps))
        rss.fit_pairwise_slopes(algorithm=args.algorithm, group_cutoff=grps)

        img_hdu = pyfits.ImageHDU(data=numpy.array(rss.cube_results[0]))
        img_hdu.header['N_GRPS'] = grps
        out_hdulist.append(img_hdu)

    pyfits.HDUList(out_hdulist).writeto(args.output_fn, overwrite=True)



