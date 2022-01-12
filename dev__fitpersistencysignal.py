#!/usr/bin/env python3

import sys
import os
import argparse
import numpy
import scipy
import scipy.optimize
import matplotlib.pyplot as plt

import astropy.io.fits as pyfits

import rss_reduce


def _persistency_plus_signal_fit_fct(p, x):
    # y = numpy.zeros(x.shape)
    # for i in range(p.shape[0]):
    #     y += p[i] * x ** (i + 1)
    y = p[0] * x + p[1] * (1. - numpy.exp(-x/p[2]))
    return y


def _persistency_plus_signal_fit_err_fct(p, x, y):
    yfit = _persistency_plus_signal_fit_fct(p, x)
    err = numpy.sqrt(y + 10 ** 2)
    return ((y - yfit) / err)


def _fit_persistency_plus_signal_pixel(_x, _y):

    good4fit = numpy.isfinite(_x) & numpy.isfinite(_y)
    x = _x[good4fit]
    y = _y[good4fit]
    if (numpy.sum(good4fit) < 5):
        # if there's no good data we can't do any fitting
        return numpy.array([1., 0., 0.])  # assume perfect linearity

    pinit = [1., 0., 1.]
    fit = scipy.optimize.leastsq(
        func=_persistency_plus_signal_fit_err_fct, x0=pinit,
        args=(x, y),
        full_output=1
    )
    # print(fit)
    pfit = fit[0]

    return pfit


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

    if (args.write_dumps):
        print("File-dumping enabled!")

    for fn in args.files:

        rss = rss_reduce.RSS(fn, max_number_files=args.max_number_files)

        if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
            rss.read_nonlinearity_corrections(args.nonlinearity_fn)
        rss.reduce(write_dumps=args.write_dumps,
                   mask_bad_data=rss_reduce.RSS.mask_SATURATED)

        # Figure out what the incremental exposure time per read is
        exptime = rss.first_header['USEREXP'] / 1000. # raw values are in milli-seconds
        delta_exptime = exptime / rss.first_header['NGROUPS']

        # now that we have the dark-rate, apply the correction to the frame to estimate the noise
        integ_exp_time = numpy.arange(rss.image_stack.shape[0]) * delta_exptime

        for (x,y) in [(1384,576), (1419,605), (1742,540), (1722,514),
            (1785,550), (1784,550), (1782,541), (1793,552), (1801,551), (1771,534), (1761,546), (1762,546),
            (1763,546), (1764,546), (1764,547), (1764,549), (1761,551), (1759,552), (1757,542), (1756,542),
            (1755,542), (1754,542), (1751,506), (1752,506), (1753,506), (1754,506),
            ]:

            _x = x-1
            _y = y-1

            bad_data = rss.bad_data_mask[:, _y, _x]

            series = rss.linearized_cube[:, _y, _x]

            # print(series.shape)

            bestfit = _fit_persistency_plus_signal_pixel(integ_exp_time[~bad_data], series[~bad_data])
            print("BESTFIT:", x,y,bestfit)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(integ_exp_time[~bad_data], series[~bad_data], s=8, c='blue')
            ax.scatter(integ_exp_time[bad_data], series[bad_data], s=4, c='grey')
            timeaxis = numpy.linspace(0, integ_exp_time[-1], 250)
            modely = _persistency_plus_signal_fit_fct(bestfit, timeaxis)
            ax.plot(timeaxis, modely)
            plot_fn = "%s____persistency_plus_signal__%04d-%04d.png" % (rss.filebase, x,y)
            fig.suptitle("F = %.0f + %.0f x (1. - exp(-x/%.3f)" % (bestfit[0], bestfit[1], bestfit[2]))
            ax.set_xlabel("integration time [seconds]")
            ax.set_ylabel("flux above read #0 [counts]")
            fig.savefig(plot_fn, dpi=200)
            plt.close(fig)

        # mean_rate_subtracted = rss.linearized_cube - integ_exp_time.reshape((-1,1,1)) * darkrate.reshape((1, darkrate.shape[0], darkrate.shape[1]))
        # print("mean rate shape:", mean_rate_subtracted.shape)
        #
        # dark_hdu = pyfits.HDUList([
        #     pyfits.PrimaryHDU(header=rss.first_header),
        #     pyfits.ImageHDU(data=darkrate, name='DARKRATE')
        # ])
        #
        # if args.output_fn is None:
        #     out_fn = rss.filebase + ".darkrate.fits"
        # else:
        #     out_fn = args.output_fn
        # print("Writing darkrate image to %s ..." % (out_fn))
        # dark_hdu.writeto(out_fn, overwrite=True)

        # for (x,y) in [(1384,576), (1419,605), (1742,540), (1722,514)]:
        #     rss.plot_pixel_curve(x,y,filebase=rss.filebase+"___")
        #     rss.dump_pixeldata(x,y,filebase=rss.filebase+"___", extras=[mean_rate_subtracted])
        #

        # rss.plot_pixel_curve(1384, 576, filebase="darkgood__" + rss.filebase+"__")
        # rss.plot_pixel_curve(1419, 605, filebase="darkgood__" + rss.filebase+"__")
        # rss.plot_pixel_curve(1742, 540, filebase="darkbad__" + rss.filebase+"__")
        # rss.plot_pixel_curve(1722, 514, filebase="darkbad__" + rss.filebase+"__")

