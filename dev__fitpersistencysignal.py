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


def _persistency_plus_signal_fit_fct(p, read_time):
    # y = numpy.zeros(x.shape)
    # for i in range(p.shape[0]):
    #     y += p[i] * x ** (i + 1)
    signal = numpy.ones_like(read_time) * p[0] + p[1] * numpy.exp(-read_time/p[2])
    return signal


def _persistency_plus_signal_fit_err_fct(p, read_time, rate, uncert):
    rate_fit = _persistency_plus_signal_fit_fct(p, read_time)
    err = uncert #numpy.sqrt(y + 10 ** 2)
    return ((rate - rate_fit) / err)


def fit_persistency_plus_signal_pixel(_x, _y, uncertainties):

    good4fit = numpy.isfinite(_x) & numpy.isfinite(_y)
    read_time = _x[good4fit]
    rate = _y[good4fit]
    uncert = uncertainties[good4fit]

    avg_rate = numpy.mean(rate)

    fallback_solution = [avg_rate, 0, 0]
    if (numpy.sum(good4fit) < 5):
        # if there's no good data we can't do any fitting
        return numpy.array(fallback_solution)  # assume perfect linearity

    # variables are: linear_rate, persistency_amplitude, persistency_timescale
    pinit = [numpy.min(rate), 2*numpy.max(rate), 3.5]
    fit = scipy.optimize.leastsq(
        func=_persistency_plus_signal_fit_err_fct, x0=pinit,
        args=(read_time, rate, uncert),
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
        # exptime = rss.first_header['USEREXP'] / 1000. # raw values are in milli-seconds
        # delta_exptime = exptime / rss.first_header['NGROUPS']

        # now that we have the dark-rate, apply the correction to the frame to estimate the noise
        # integ_exp_time = numpy.arange(rss.image_stack.shape[0]) * delta_exptime

        for (x,y) in [(1384,576), (1419,605), (1742,540), (1722,514),
            (1785,550), (1784,550), (1782,541), (1793,552), (1801,551), (1771,534), (1761,546), (1762,546),
            (1763,546), (1764,546), (1764,547), (1764,549), (1761,551), (1759,552), (1757,542), (1756,542),
            (1755,542), (1754,542), (1751,506), (1752,506), (1753,506), (1754,506),
            ]:

            _x = x-1
            _y = y-1

            bad_data = rss.bad_data_mask[:, _y, _x]
            # series = rss.linearized_cube[:, _y, _x]
            rate_series = rss.differential_stack[:, _y, _x]

            # TODO: implement better noise model, accounting for read-noise and gain
            uncertainties = numpy.sqrt(rss.image_stack[:, _y, _x])
            # print(series.shape)

            numpy.savetxt(
                "persistency_dump_%04dx%04d.txt" % (x,y),
                numpy.array([rss.read_times, rate_series, uncertainties, bad_data]).T,
            )

            bestfit = fit_persistency_plus_signal_pixel(
                rss.read_times[~bad_data], rate_series[~bad_data], uncertainties[~bad_data]
            )
            print("BESTFIT:", x,y,bestfit)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(rss.read_times[~bad_data], rate_series[~bad_data], s=8, c='blue')
            ax.scatter(rss.read_times[bad_data], rate_series[bad_data], s=4, c='grey')
            timeaxis = numpy.linspace(0, numpy.nanmax(rss.read_times), 250)
            # print("read-times = \n", timeaxis)
            modely = _persistency_plus_signal_fit_fct(bestfit, timeaxis)
            # print("best-fit:", bestfit)
            # print("model-fit = \n", modely)

            ax.plot(timeaxis, modely)
            plot_fn = "%s____persistency_plus_signal__%04d-%04d.png" % (rss.filebase, x,y)
            fig.suptitle("F = %.0f + %.0f x exp(-t/%.3fs)" % (bestfit[0], bestfit[1], bestfit[2]))
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

