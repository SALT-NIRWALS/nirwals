#!/usr/bin/env python3

import sys
import os
import numpy
import astropy.io.fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
import argparse

import nirwals_reduce


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
    cmdline.add_argument("--coords", nargs="+",
                         help="list of pixel coordinates")
    args = cmdline.parse_args()

    if (args.write_dumps):
        print("File-dumping enabled!")

    for fn in args.files:

        rss = rss_reduce.RSS(fn)
        # rss.reduce(write_dumps=False)
        # rss.write_results()

        x = int(sys.argv[2])
        y = int(sys.argv[3])
        _x = x - 1
        _y = y - 1

        print("Starting image input")
        rss.load_all_files()
        # rss.subtract_first_read()

        pyfits.PrimaryHDU(data=rss.image_stack[35]).writeto("dump_frame.fits", overwrite=True)

        print("Extracting data")

        # remove zero
        raw_series = numpy.array(rss.image_stack[:, _y,_x])
        # rss.subtract_first_read()

        # extract data
        series = raw_series - raw_series[0]

        diffs = numpy.pad(numpy.diff(raw_series), (1,0))
        print(diffs.shape, series.shape)

        # flag bad/saturated pixels
        max_diff = numpy.nanpercentile(diffs, 90)
        bad = (raw_series > 63000) | (diffs < 0.3 * max_diff)

        # numpy.savetxt("pixeldump_x%04d_y%04d.raw" % (x,y), raw_series)
        # numpy.savetxt("pixeldump_x%04d_y%04d.txt" % (x,y), series)

        avg_rate = series[10]/10.
        print(avg_rate)

        integrations_count = numpy.arange(series.shape[0])

        pfit2 = rss._fit_nonlinearity_pixel(integrations_count[~bad], series[~bad])
        print(pfit2)

        best_fit_direct = pfit2[0] * integrations_count + \
                          pfit2[1] * numpy.power(integrations_count,2) + \
                          pfit2[2] * numpy.power(integrations_count,3)

        integrations_count = numpy.arange(series.shape[0])
        computed_countrate = integrations_count * pfit2[0]
        # 1.00000000e+00 2.36453272e-07 1.13436118e-10] = avg_rate

        print("Fitting pixel")
        pfit = rss._fit_nonlinearity_pixel(series[~bad], computed_countrate[~bad])
        print(pfit)

        linearized = pfit[0] * numpy.power(series, 1) + \
            pfit[1] * numpy.power(series, 2) + \
            pfit[2] * numpy.power(series, 3)

        numpy.savetxt("pixeldump_x%04d_y%04d.complete" % (x,y),
            numpy.array([integrations_count, raw_series,
                         series, computed_countrate, best_fit_direct, linearized,
                         diffs, bad.astype(numpy.int)]).T
        )

        last_good_sample = numpy.max(integrations_count[~bad])
        print("last good sample:", last_good_sample)
        # make plot with all information
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.scatter(integrations_count, raw_series, s=4, label="raw data")
        ax.scatter(integrations_count, series, s=8, label="zero subtracted")
        ax.scatter(integrations_count[bad], series[bad], c='grey', marker='x', s=16, label='bad', linewidth=1)
        ax.plot(integrations_count, best_fit_direct, label="fit")
        ax.plot(integrations_count, computed_countrate, label='constant rate')
        ax.scatter(integrations_count, linearized, s=8, label='non-lin corrected')
        ax.scatter(integrations_count, linearized + raw_series[0], s=3)
        ax.legend()
        ax.set_ylim((-500, 74000))
        ax.set_xlim((-0.5, numpy.max(integrations_count) + 2.5))
        ax.axhline(y=63000, linestyle=':', color='grey')
        ax.axvline(x=last_good_sample + 0.5, linestyle=":", color='grey')
        ax.set_xlabel("Read")
        ax.set_ylabel("counts [raw/corrected]")
        fig.suptitle("x=%d  y=%d" % (x,y))
        fig.tight_layout()
        fig.show()
        fig.savefig("single_pixel__x%04d_y%04d.png" % (x,y), dpi=300)

        # rss.fit_nonlinearity(ref_frame_id=4)

        # rss.plot_pixel_curve(818,1033)
        # rss.plot_pixel_curve(1700,555)
