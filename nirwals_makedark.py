#!/usr/bin/env python3

import sys
import os
import argparse
import numpy

import astropy.io.fits as pyfits

import nirwals_reduce

from dev__fitpersistencysignal import fit_persistency_plus_signal_pixel

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

        # count the number of good samples
        number_good_dark_samples = numpy.sum(~rss.bad_data_mask, axis=0)
        print("# good samples", number_good_dark_samples.shape)

        # Figure out what the incremental exposure time per read is
        exptime = rss.first_header['USEREXP'] / 1000. # raw values are in milli-seconds
        delta_exptime = exptime / rss.first_header['NGROUPS']
        integ_exp_time = numpy.arange(rss.image_stack.shape[0]) * delta_exptime

        darkrate = rss.weighted_mean / delta_exptime

        # keep track of what kind of pixel each one is
        pixeltype = numpy.full_like(darkrate, fill_value=rss_reduce.darktype_GOOD, dtype=numpy.int)

        # find typical noise levels, so we can identify pixels that may need special treatment
        good_data = numpy.isfinite(darkrate)
        for iter in range(3):
            _stats = numpy.nanpercentile(darkrate[good_data], [16, 50, 84])
            print(iter, _stats)
            _median = _stats[1]
            _sigma = 0.5 * (_stats[2]-_stats[0])
            good_data = good_data & (darkrate > _median-3*_sigma) & (darkrate < _median+3*_sigma)
        final_stats = numpy.nanpercentile(darkrate[good_data], [16, 50, 84])
        print("final image stats: ", final_stats)

        darkrate_cleaned = numpy.copy(darkrate)
        # darkrate_cleaned[~good_data] = numpy.NaN
        pyfits.PrimaryHDU(data=darkrate_cleaned).writeto("darkrate_cleaned.fits", overwrite=True)

        # identify pixels we need to have a special treatment for

        negative_pixels = ~good_data & (darkrate < 0)
        # let's ignore those for now, and set the actual dark-current to 0
        # TODO: add code here
        pixeltype[negative_pixels] = rss_reduce.darktype_COLD

        min_exptime_required = 15. # seconds
        min_frames = min_exptime_required / (delta_exptime)
        hot_pixels = ~good_data & (darkrate > 0) & (number_good_dark_samples < min_frames)
        pixeltype[hot_pixels] = rss_reduce.darktype_HOT

        warm_pixels = ~good_data & (darkrate > 0) & (number_good_dark_samples >= min_frames)
        pixeltype[warm_pixels] = rss_reduce.darktype_WARM
        iy,ix = numpy.indices(darkrate.shape)
        ixy = numpy.dstack([ix,iy])
        print("IXY", ixy.shape)
        warm_ixy = ixy[warm_pixels]
        print("warm_ixy:", warm_ixy.shape)

        warm_coefficients = numpy.zeros((3, darkrate.shape[0], darkrate.shape[1]))
        for i, _xy in enumerate(warm_ixy):
            x, y = _xy[0], _xy[1]
            fullseries = rss.linearized_cube[:, y, x]
            bad = rss.bad_data_mask[:, y, x]

            series = fullseries[~bad]
            times = integ_exp_time[~bad]
            bestfit = fit_persistency_plus_signal_pixel(times, series)

            warm_coefficients[:, y,x] = bestfit
            if ((i%1000) == 0):
                sys.stdout.write(">")
                sys.stdout.flush()

        # now that we have the dark-rate, apply the correction to the frame to estimate the noise
        mean_rate_subtracted = rss.linearized_cube - integ_exp_time.reshape((-1,1,1)) * darkrate.reshape((1, darkrate.shape[0], darkrate.shape[1]))
        print("mean rate shape:", mean_rate_subtracted.shape)

        dark_hdu = pyfits.HDUList([
            pyfits.PrimaryHDU(header=rss.first_header),
            pyfits.ImageHDU(data=darkrate, name='DARKRATE'),
            pyfits.ImageHDU(data=number_good_dark_samples, name='N_SAMPLES'),
            pyfits.ImageHDU(data=warm_coefficients, name='WARM_DARK_COEFFS'),
        ])

        if args.output_fn is None:
            out_fn = rss.filebase + ".darkrate.fits"
        else:
            out_fn = args.output_fn
        print("Writing darkrate image to %s ..." % (out_fn))
        dark_hdu.writeto(out_fn, overwrite=True)

        for (x,y) in [(1384,576), (1419,605), (1742,540), (1722,514)]:
            rss.plot_pixel_curve(x,y,filebase=rss.filebase+"___")
            rss.dump_pixeldata(x,y,filebase=rss.filebase+"___", extras=[mean_rate_subtracted])


        # rss.plot_pixel_curve(1384, 576, filebase="darkgood__" + rss.filebase+"__")
        # rss.plot_pixel_curve(1419, 605, filebase="darkgood__" + rss.filebase+"__")
        # rss.plot_pixel_curve(1742, 540, filebase="darkbad__" + rss.filebase+"__")
        # rss.plot_pixel_curve(1722, 514, filebase="darkbad__" + rss.filebase+"__")

