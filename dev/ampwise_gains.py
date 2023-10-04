#!/usr/bin/env python3

import os
import astropy.io.fits as pyfits
import scipy
import numpy
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.stats

import multiparlog as mplog
import logging
import argparse
from matplotlib.backends.backend_pdf import PdfPages

def listhandler(inlist):

    logger = logging.getLogger('ListHandler')
    outlist = []
    for fn in inlist:
        if (fn.startswith("@") and os.path.isfile(fn[1:])):
            addlist = []
            with open(fn[1:], "r") as f:
                lines = f.readlines()
                for l in lines:
                    if l.strip().startswith("#"):
                        continue
                    addlist.append(l.strip().split()[0])
            outlist.extend(listhandler(addlist))
        elif os.path.isfile(fn):
            logger.info("Adding %s to input filelist" % (fn))
            outlist.append(fn)
        else:
            logger.debug("Input (%s) not found or understood" % (fn))

    return outlist

if __name__ == "__main__":


    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("Gain")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--y1", dest="y1", default=0, type=int,
                         help="start row to use for stacking")
    cmdline.add_argument("--y2", dest="y2", default=0, type=int,
                         help="end row to use for stacking")
    cmdline.add_argument("--master", dest="master", default=False, action='store_true',
                         help="combine results for all input frames")
    cmdline.add_argument("files", nargs='+', help="reduced input files")
    args = cmdline.parse_args()

    logger.info("Working dir: %s" % (os.getcwd()))

    raw_filelist = args.files
    filelist = listhandler(raw_filelist)

    master_ratios = []
    master_gains = []
    for fn in filelist:

        logger.info("Loading input file (%s)" % (fn))
        hdulist = pyfits.open(fn)
        data = hdulist['SCI'].data
        # print(data.shape)

        pdf_fn = fn[:-5]+"__gainfit.pdf"
        pdf = PdfPages("gainfit2.pdf")

        namps = hdulist[0].header['NOUTPUTS']
        logger.info("Preparing for %d amplifiers" % (namps))

        logger.info("Stacking rows (%d ... %d" % (args.y1, args.y2))
        colstack = numpy.sum(data[args.y1:args.y2], axis=0)

        # find only good data -- first & last 4 colums are reference pixels and always bad
        good_data = numpy.isfinite(colstack)
        good_data[:4] = False
        good_data[-4:] = False

        namps = hdulist[0].header['NOUTPUTS']
        ampsize = 2048 // namps

        fig,ax = plt.subplots()
        ax.scatter(numpy.arange(colstack.shape[0]), colstack, s=1)
        for n in range(namps+1):
            ax.axvline(x=n*ampsize, alpha=0.2)
        pdf.savefig(fig)

        amp_x = numpy.arange(ampsize)
        n_points = 10
        point_spacing = ampsize // n_points

        chunks = []
        amp_edges = []
        for amp in range(namps):

            a1 = amp * ampsize
            a2 = a1 + ampsize

            amp_good = good_data[a1:a2]
            amp_x = numpy.arange(ampsize)[amp_good]
            amp_val = colstack[a1:a2][amp_good]

            point_spacing = (numpy.max(amp_x) - numpy.min(amp_x)) / n_points
            base_points = (numpy.arange(n_points, dtype=float) + 0.5) * point_spacing + numpy.min(amp_x)

            # do fit
            lsq = scipy.interpolate.LSQUnivariateSpline(
                x=amp_x, y=amp_val,
                t=base_points,
                k=3,
                ext='extrapolate')

            linfit = scipy.stats.linregress(x=amp_x, y=amp_val)

            amp_edges.append([linfit.intercept, linfit.intercept+linfit.slope*ampsize])

            # print(amp, lsq(0), lsq(ampsize))
            logger.debug("%02d %10.1f %10.1f" % (amp, linfit.intercept, linfit.intercept+linfit.slope*ampsize))
            chunks.append( (amp_x+a1, linfit.slope*amp_x+linfit.intercept) )

        fig,ax = plt.subplots(figsize=(13,5))
        ax.scatter(numpy.arange(colstack.shape[0]), colstack, s=1)
        for n in range(namps+1):
            ax.axvline(x=n*ampsize, alpha=0.2)
        for [x,y] in chunks:
            ax.plot(x,y)
        pdf.savefig(fig)


        amp_edges = numpy.array(amp_edges)
        ratios = amp_edges[1:,0] / amp_edges[:-1,1]
        logger.debug(ratios)


        gains = [1.]
        _gain = 1.0
        for g in ratios:
            _gain = _gain * g
            gains.append(_gain)

        gains = numpy.array(gains)
        logger.debug(gains.shape)


        gains_full = numpy.repeat(gains, ampsize)
        logger.debug(gains_full.shape)

        dump_fn = fn[:-5]+"__gaindump.txt"
        numpy.savetxt(dump_fn, gains.reshape((-1,1)))

        img_fixed = data / gains_full.reshape((1,-1))
        out_fn = fn[:-5]+"__gainfixed.fits"
        logger.info("Writing output to %s" % (out_fn))
        pyfits.PrimaryHDU(data=img_fixed).writeto(out_fn, overwrite=True)

        # keep track of all ratios & gains
        master_gains.append(gains)
        master_ratios.append(ratios)


    master_gains = numpy.array(master_gains)
    master_ratios = numpy.array(master_ratios)

    pyfits.PrimaryHDU(data=master_gains).writeto("mastergains.fits", overwrite=True)
    pyfits.PrimaryHDU(data=master_ratios).writeto("masterratios.fits", overwrite=True)
