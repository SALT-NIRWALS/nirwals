#!/usr/bin/env python3

import astropy.io.fits as pyfits
import scipy
import numpy
import matplotlib.pyplot as plt
import scipy.interpolate
import multiparlog as mplog
import logging
import argparse
from matplotlib.backends.backend_pdf import PdfPages

if __name__ == "__main__":


    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("Gain")

    cmdline = argparse.ArgumentParser()
    # cmdline.add_argument("--stage", dest="staging_dir", default="./",
    #                      help="staging directory")
    # cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
    #                      help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--y1", dest="y1", default=0, type=int,
                         help="start row to use for stacking")
    cmdline.add_argument("--y2", dest="y2", default=0, type=int,
                         help="end row to use for stacking")
    # cmdline.add_argument("--test", dest="test", default=None,
    #                      help="test mode; syntax: --test=delay:@filelist")
    cmdline.add_argument("files", nargs=1, help="reduced input files")
    args = cmdline.parse_args()

    fn = args.files[0]
    logger.info("Loading input file (%s)" % (fn))
    hdulist = pyfits.open(fn)
    data = hdulist['SCI'].data
    # print(data.shape)

    namps = hdulist[0].header['NOUTPUTS']
    logger.info("Preparing for %d amplifiers" % (namps))

    logger.info("Stacking rows (%d ... %d" % (args.y1, args.y2))
    colstack = numpy.sum(data[args.y1:args.y2], axis=0)

    # find only good data -- first & last 4 colums are reference pixels and always bad
    good_data = numpy.isfinite(colstack)
    good_data[:4] = False
    good_data[-4:] = False

    xvals = numpy.arange(colstack.shape[0])
    colstack4fit = colstack[good_data]
    xvals4fit = xvals[good_data]

    n_basepoints=[4,4,8,16,12]
    n_iterations = len(n_basepoints)
    ampsize = 2048 // namps

    pdf = PdfPages("gainfit.pdf")

    for iteration in range(n_iterations):

        logger.info("Starting on iteration %d of %d" % (iteration+1, n_iterations))
        # find basepoints
        n_points = n_basepoints[iteration]
        point_spacing = 2048 / n_points
        base_points = (numpy.arange(n_points, dtype=float) + 0.5) * point_spacing

        # do fit
        lsq = scipy.interpolate.LSQUnivariateSpline(
            x=xvals4fit, y=colstack4fit,
            t=base_points,
            k=3,
            ext='zeros')

        raw_correction = colstack/lsq(xvals)
        corr = numpy.nanmean(raw_correction.reshape((namps,-1)), axis=1)

        corr_full = numpy.repeat(corr, ampsize)

        corr_colstack = colstack / corr_full

        fig,axs = plt.subplots(ncols=2, figsize=(10,4))
        fig.suptitle("iteration %d" % (iteration+1))
        xval = numpy.arange(colstack.shape[0])
        ax = axs[0]
        ax.scatter(xvals4fit, colstack4fit, s=1)
        ax.scatter(xvals[good_data], corr_colstack[good_data], s=1)
        ax.plot(xvals[good_data], lsq(xvals)[good_data], color='orange', alpha=0.7)
        ax.scatter(base_points, lsq(base_points), s=40, alpha=1, color='grey')

        axs[1].scatter(xvals, colstack/lsq(xvals), s=0.3)
        axs[1].set_ylim((0.8,1.2))
        colstack4fit = corr_colstack[good_data]

        pdf.savefig(fig)

    scaling_factors = lsq(xval) / colstack
    corr = numpy.nanmean(scaling_factors.reshape((namps,-1)), axis=1)
    fig,ax = plt.subplots(figsize=(10,3))
    fig.suptitle("final gain corrections")
    ax.scatter(xval, scaling_factors, s=0.4)
    ax.scatter( (numpy.arange(corr.shape[0], dtype=float)+0.5)*ampsize, corr )
    for n in range(namps):
        ax.axvline(x=n*ampsize, alpha=0.2, color='grey')
    print(corr)
    pdf.savefig(fig)

    pdf.close()

    logger.info("all done!")
