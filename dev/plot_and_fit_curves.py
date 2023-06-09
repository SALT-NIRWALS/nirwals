#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy
import pandas
import glob
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
import scipy.optimize
import sys
import os

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

def _persistency_plus_signal_fit_fct2(p, read_time):
    # y = numpy.zeros(x.shape)
    # for i in range(p.shape[0]):
    #     y += p[i] * x ** (i + 1)
    signal = numpy.ones_like(read_time) * p[0] + p[1] * numpy.exp(-read_time/p[2]) + p[3] * numpy.exp(-read_time/p[4])
    return signal
def _persistency_plus_signal_fit_err_fct2(p, read_time, rate, uncert):
    rate_fit = _persistency_plus_signal_fit_fct2(p, read_time)
    err = uncert #numpy.sqrt(y + 10 ** 2)
    return ((rate - rate_fit) / err)

if __name__ == "__main__":

    filename = sys.argv[1]
    pixel_fn = sys.argv[2]

    dir, _bn = os.path.split(filename)
    bn,ext = os.path.splitext(_bn)

    pdf_fn = "%s.pdf" % (bn)
    fit_txt_fn = "%s_fitresults.txt" % (bn)

    with open(pixel_fn, "r") as pf:
        lines = pf.readlines()
        xy = [[int(i) for i in line.split()[:2]] for line in lines]
    xy = numpy.array(xy)
    # print(xy)


    data = pandas.read_csv(filename, nrows=1, delim_whitespace=True)
    n_cols = len(data.columns)-1
    col_names = ['pixel_%02d' % (i+1) for i in range(n_cols)]

    print("Reading data from %s" % (filename))
    data = pandas.read_csv(filename, delim_whitespace=True, names=['read'] + col_names, index_col='read')
#     data.info()

    print("Saving plots to %s" % (pdf_fn))
    pdf = PdfPages(pdf_fn)
    final_fit_results = []

    for icol, colname in enumerate(col_names[:]):
        # print(colname)
        pixel_xy = xy[icol]
        pixel_xy_str = "pixel @ %d , %d" % (pixel_xy[0], pixel_xy[1])

        sys.stdout.write("\rWorking on pixel % 4d of % 4d (%04d / %04d)" % (icol+1, len(col_names), pixel_xy[0], pixel_xy[1]))
        sys.stdout.flush()

        dt = numpy.diff(data.index)
        flux = numpy.array(data[colname])
        rate = flux[1:] - flux[:-1]
        rate2 = numpy.diff(data[colname])
        uncert = 25.

        valid_data = flux[1:] < 33000

        rate_guess = 0. #numpy.max([0., numpy.min([60000, numpy.nanmin(rate)])])
        pers_guess = 100. #numpy.max([0, numpy.min([60e3, (numpy.nanmax(rate) - rate_guess)])])
        # tau: typical value found in data fits
        tau_guess = 100.
        pinit = [rate_guess, pers_guess, tau_guess]
        # fit_results = scipy.optimize.least_squares(
        #     fun=_persistency_plus_signal_fit_err_fct,
        #     x0=pinit,
        #     bounds=([0, 0, 1.], [numpy.Inf, 65e3, 1000]),
        #     kwargs=dict(read_time=data.index[1:][valid_data], rate=rate[valid_data],
        #                 uncert=numpy.full_like(rate, fill_value=30)[valid_data]),
        # )

        pinit = [rate_guess, pers_guess, 2.0, 0, 50.]
        fit_results = scipy.optimize.least_squares(
            fun=_persistency_plus_signal_fit_err_fct2,
            x0=pinit,
            bounds=([0, 0, 1., 0, 5], [numpy.Inf, 65e3, 4, 65e3, 100]),
            kwargs=dict(read_time=data.index[1:][valid_data], rate=rate[valid_data],
                        uncert=numpy.full_like(rate, fill_value=30)[valid_data]),
        )


        bestfit = fit_results.x
        bounds_limited_mask = fit_results.active_mask
        fit_successful = fit_results.success

        rate_fit = _persistency_plus_signal_fit_fct2(p=bestfit, read_time=data.index[1:])
        cumulative_rate_fit = numpy.cumsum(rate_fit)

        pers_only_fit = _persistency_plus_signal_fit_fct2(p=[0, bestfit[1], bestfit[2], 0., 0.],
                                                         read_time=data.index[1:])
        signal_only_fit = _persistency_plus_signal_fit_fct2(p=[bestfit[0], 0., 0., 0., 0.],
                                                         read_time=data.index[1:])
        longpers_only_fit = _persistency_plus_signal_fit_fct2(p=[0, 0, 0, bestfit[3], bestfit[4]],
                                                              read_time=data.index[1:])
        full_persistency = _persistency_plus_signal_fit_fct2(p=[0, bestfit[1], bestfit[2], bestfit[3], bestfit[4]],
                                                              read_time=data.index[1:])
        pers_only_cumulative_rate_fit = numpy.cumsum(pers_only_fit)
        signal_only_cumulative_rate_fit = numpy.cumsum(signal_only_fit)
        longpers_only_cumulative_rate_fit = numpy.cumsum(longpers_only_fit)
        fullpers_only_cumulative_rate_fit = numpy.cumsum(full_persistency)

        # print(data.index.shape,
        #     pers_only_fit.shape, signal_only_fit.shape,
        #         pers_only_cumulative_rate_fit.shape,
        #       signal_only_cumulative_rate_fit.shape,
        #      data[colname].shape
        #      )

        fig = plt.figure(figsize=(8,4))

        ax = fig.add_subplot(121)
        ax.set_title("%s \n signal=%.4f" % (pixel_xy_str, bestfit[0]))
        ax.scatter(data.index, data[colname], s=1, c='blue', label='raw')
        ax.plot(data.index[1:], cumulative_rate_fit, c='blue') #label='fit')

        # print(data.index[1:].shape, flux[1:].shape, pers_only_cumulative_rate_fit.shape)
        # print(data.shape)
        ax.scatter(data.index[1:], flux[1:]-fullpers_only_cumulative_rate_fit,
                       s=0.5, c='orange')
        # ax.scatter(data.index[1:], flux[1:]-signal_only_cumulative_rate_fit,
        #                s=0.5, c='green',)

        ax.plot(data.index[1:], signal_only_cumulative_rate_fit, c='orange', label='signal')
        ax.plot(data.index[1:], pers_only_cumulative_rate_fit, c='green', label='short pers.')
        ax.plot(data.index[1:], longpers_only_cumulative_rate_fit, c='red', label='long pers.')

        ax.legend(loc='lower right', markerscale=5., fontsize=6)
        [_imin,_imax] = ax.get_ylim()
        _min = numpy.max([-100, _imin])
        _max = numpy.min([_imax, 45000])
        ax.set_ylim([_min, _max])


        ax2 = fig.add_subplot(122)
        ax2.scatter(data.index[1:], rate, s=2)
        # ax2.scatter(data.index[1:], rate2, s=2)
        ax2.plot(data.index[1:], rate_fit)
        ax2.plot(data.index[1:], pers_only_fit, color='green')
        ax2.plot(data.index[1:], longpers_only_fit, color='red')

        ax2.axhline(y=bestfit[0], c='orange')
        ax2.set_title("short: tau=%.2fs  amp=%.1f cts/s\nfast: tau=%.1f  amp=%.1f" % (bestfit[2], bestfit[1], bestfit[4], bestfit[3]))

        pdf.savefig(fig)
        plt.close(fig)

        final_fit_results.append([pixel_xy[0], pixel_xy[1], bestfit[0], bestfit[1], bestfit[2], bestfit[3], bestfit[4]])

    pdf.close()

    print("\n")
    print("Saving fit results to %s" % (fit_txt_fn))
    numpy.savetxt(fit_txt_fn, numpy.array(final_fit_results))

    print("all done!")
