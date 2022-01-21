#!/usr/bin/env python3
import multiprocessing
import os
import sys
import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import itertools
import multiprocessing
import argparse

from astropy import log
log.setLevel('ERROR')

import astropy
print(astropy.__path__)

# def fit_nonlinearity_sequence(pinit, args):
#     """
#     Mostly a wrapper around scipy.optimize.leastsq.
#
#     """
#
#     def fit_fct(p, x):
# #         y = numpy.zeros(x.shape)
# #         for i in range(p.shape[0]):
# #             y += p[i] * x**(i+1)
# #         return y
#     def err_fct(p,x,y,err, fitrange_x, fitrange_y):
#         yfit = fit_fct(p,x)
#         in_fit_range = numpy.isfinite(x) & numpy.isfinite(y)
#         if (fitrange_x is not None):
#             in_fit_range = in_fit_range & (x >= fitrange_x[0]) & (x < fitrange_x[1])
#         if (fitrange_y is not None):
#             in_fit_range = in_fit_range & (y >= fitrange_y[0]) & (y < fitrange_y[1])
#         if (err is None):
#             return ((y-yfit))[in_fit_range]
#         return ((y-yfit)/err)[in_fit_range]
#
#     # (medlevel, exptime, None, intensity_range, exptime_range) = args
#
#     fit = scipy.optimize.leastsq(err_fct, pinit, args=args, full_output=1)
#
#     pfit = fit[0]
#
#     uncert = numpy.sqrt(numpy.diag(fit[1]))
#     print(pfit, uncert)
#     return pfit, uncert
#

darktype_GOOD = 0
darktype_COLD = 1
darktype_WARM = 2
darktype_HOT = 3

class RSS(object):

    mask_SATURATED = 0x0001
    mask_LOW_RATE = 0x0002
    mask_BAD_DARK = 0x0004
    mask_NEGATIVE = 0x0008

    def __init__(self, fn, max_number_files=-1,
                 saturation_level=60000,
                 saturation_fraction=0.25, saturation_percentile=95):
        self.fn = fn
        self.filelist = []

        self.image_stack_initialized = False
        self.first_read_subtracted = False
        self.first_read = None
        self.first_header = None

        self.nonlin_fn = None
        self.nonlinearity_cube = None

        # store values we may/will need during reduction
        self.max_number_files = -1 if max_number_files is None else max_number_files
        self.saturation_level = saturation_level
        self.saturation_fraction = saturation_fraction
        self.saturation_percentile = saturation_percentile

        self.get_full_filelist()

    def get_full_filelist(self):
        # get basedir
        fullpath = os.path.abspath(self.fn)
        print(fullpath)
        self.basedir, _fn = os.path.split(fullpath)
        self.filebase = ".".join(_fn.split(".")[:-2])
        # print(self.basedir, filebase)

        for _read in range(1,100):
            filename = "%s.%d.fits" % (self.filebase, _read)
            full_filename = os.path.join(self.basedir, filename)
            if (os.path.isfile(full_filename)):
                self.filelist.append(full_filename)
            else:
                break

            # print(full_filename, os.path.isfile(full_filename))

        print("\n -- ".join(self.filelist))

        return

    def add_file(self, filename):
        return

    def store_header_info(self, hdr):
        self.first_header = hdr

        self.exptime = hdr['USEREXP'] / 1000.
        self.n_groups = hdr['NGROUPS']
        self.diff_exptime = self.exptime / self.n_groups

    def load_all_files(self, max_number_files=None):

        if (max_number_files is None):
            max_number_files = self.max_number_files
        if (self.image_stack_initialized):
            print("stack already initialized, skipping repeat try")
            return

        self._image_stack = []

        # open all frames
        _filelist = self.filelist
        if (max_number_files > 0):
            print("Limiting filelist to %d frames" % (max_number_files))
            _filelist = _filelist[:max_number_files]
        for fn in _filelist:
            hdulist = pyfits.open(fn)
            # hdulist.info()
            imgdata = hdulist[0].data
            self._image_stack.append(imgdata)
            if  (self.first_header is None):
                self.store_header_info(hdulist[0].header)

            # break

        # calculate the initial image stack
        self.image_stack = numpy.array(self._image_stack, dtype=numpy.float32)
        print(self.image_stack.shape)

        self.image_stack_initialized = True

    def reduce(self, dark_fn=None, write_dumps=False, mask_bad_data=None):

        self.load_all_files()

        # self.subtract_first_read()
        # apply first-read subtraction
        self.reset_frame = self.image_stack[0]
        reset_frame_subtracted = self.image_stack - self.reset_frame

        # apply any necessary corrections for nonlinearity and other things
        print("Applying non-linearity corrections")
        linearized = self.apply_nonlinearity_corrections(reset_frame_subtracted)
        # print("linearized = ", linearized)
        if (linearized is None):
            print("No linearized data found, using raw data instead")
            linearized = reset_frame_subtracted
        self.linearized_cube = linearized

        dark_cube = numpy.zeros_like(linearized)
        if (dark_fn is None):
            print("No dark correction requested, skipping")
        elif (not os.path.isfile(dark_fn)):
            print("Dark requested (%s) but not found" % (dark_fn))
        else:
            try:
                # load dark-rate image
                print("Loading dark-corrections from %s" % (dark_fn))
                dark_hdu = pyfits.open(dark_fn)
                dark = dark_hdu['DARKRATE'].data

                # perform a dark subtraction;
                # dark-current = rate [in cts/sec] * frame-# * exposure-time per frame [in sec]
                dark_cube = (numpy.arange(linearized.shape[0], dtype=numpy.float).reshape((-1,1,1)) + 1) \
                            * self.diff_exptime \
                            * dark.reshape((1, dark.shape[0], dark.shape[1]))
                print("shape of dark cube: ", dark_cube.shape)
                linearized -= dark_cube
            except Exception as e:
                print("error during dark subtraction:\n",e)

        # calculate differential stack
        self.differential_stack = numpy.pad(numpy.diff(linearized, axis=0), ((1,0),(0,0),(0,0)))
        print("diff stack:", self.differential_stack.shape)

        # mask out all saturated and/or otherwise bad samples
        max_count_rates = -1000 # TODO: FIX THIS numpy.nanpercentile(self.differential_stack, q=self.saturation_percentile, axis=0)
        # print("max counrates:", max_count_rates.shape)

        # TODO: implement full iterative outlier rejection here
        print("Identifying bad pixels")
        bad_data = numpy.zeros_like(self.image_stack, dtype=numpy.bool)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_SATURATED) > 0):
            bad_data = bad_data | (self.image_stack > self.saturation_level)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_LOW_RATE) > 0):
            bad_data = bad_data | (self.differential_stack < self.saturation_fraction*max_count_rates)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_BAD_DARK) > 0):
            bad_data = bad_data | (dark_cube >= linearized)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_NEGATIVE) > 0):
            bad_data = bad_data | (linearized < 0)

        self.bad_data_mask = bad_data

        # bad_data = (self.image_stack > self.saturation_level) | \
        #            (self.differential_stack < self.saturation_fraction*max_count_rates) | \
        #            (dark_cube >= linearized) | \
        #            (linearized < 0)

        print("Cleaning image cube")
        self.clean_stack = self.differential_stack.copy()
        self.clean_stack[bad_data] = numpy.NaN
        self.clean_stack[0, :, :] = numpy.NaN # mask out the first slice, which is just padding

        # calculate a average countrate image
        print("calculating final image from stack")
        # image7 = numpy.nanmean(self.clean_stack[:7], axis=0)
        self.reduced_image_plain = numpy.nanmean(self.clean_stack, axis=0)
        noise = numpy.sqrt(self.image_stack)
        noise[bad_data] = numpy.NaN
        noise[0, :, :] = numpy.NaN
        self.inv_noise = numpy.nansum(1./noise, axis=0)
        self.weighted_mean = numpy.nansum(self.clean_stack / noise, axis=0) / self.inv_noise
        self.noise_image = 1. / self.inv_noise
        # print(image.shape)

        # ratios = linearized / linearized[3:4, :, :]

        if (write_dumps):
            print("Writing all dumps")
            bn = self.filebase + "__"
            print("Dump-file basename: ", bn)
            pyfits.PrimaryHDU(data=self.image_stack).writeto(bn+"stack_raw.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.image_stack).writeto(bn+"stack_zerosub.fits", overwrite=True)
            pyfits.PrimaryHDU(data=linearized).writeto(bn+"stack_linearized.fits", overwrite=True)
            print("writing darkcube")
            pyfits.PrimaryHDU(data=dark_cube).writeto(bn+"stack_darkcube.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.differential_stack).writeto(bn+"stack_diff.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.clean_stack).writeto(bn+"stack_clean.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=ratios).writeto("stack_ratios.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=self.reduced_image_plain).writeto("final_image.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=image7).writeto("final_image7.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=max_count_rates).writeto("max_count_rates.fits", overwrite=True)
            pyfits.PrimaryHDU(data=bad_data.astype(numpy.int)).writeto(bn+"bad_data.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.weighted_mean).writeto(bn+"final_image_weighted.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.inv_noise).writeto(bn+"final_inv_noise.fits", overwrite=True)

        return

    def subtract_first_read(self):
        if (self.first_read_subtracted):
            print("First read already subtracted, skipping")
            return

        if (self.first_read is None):
            self.first_read = self.image_stack[0].copy()

        self.image_stack -= self.first_read
        self.first_read_subtracted = True

    def _nonlinearity_fit_fct(self, p, x):
        y = numpy.zeros(x.shape)
        for i in range(p.shape[0]):
            y += p[i] * x**(i+1)
        return y

    def _nonlinearity_fit_err_fct(self, p, x, y):
        yfit = self._nonlinearity_fit_fct(p, x)
        err = numpy.sqrt(y + 10 ** 2)
        return ((y - yfit) / err)

    def _fit_nonlinearity_pixel(self, _x, _y):
        # print("fitting", _x.shape, _y.shape)
        # return 1

        # define two sub-routines we need for the fitting


        # prepare the data and initial fitting parameters
        # print(ix, iy)
        # _x = masked[:, iy, ix]
        # _y = linearized_intensity[:, iy, ix]
        good4fit = numpy.isfinite(_x) & numpy.isfinite(_y)
        # print(_x)
        # print(_y)
        x = _x[good4fit]
        y = _y[good4fit]
        if (numpy.sum(good4fit) < 5):
            # if there's no good data we can't do any fitting
            return numpy.array([1., 0., 0.]) # assume perfect linearity

        pinit = [1., 0., 0.]
        fit = scipy.optimize.leastsq(
            func=self._nonlinearity_fit_err_fct, x0=pinit,
            args=(x, y),
            full_output=1
        )
        pfit = fit[0]

        return pfit

    def fit_nonlinearity(self, ref_frame_id=10, max_linear=50000, make_plot=False):

        # self.subtract_first_read()
        # if (self.first_read_subtracted):
        #     bad_data = (self.image_stack + self.first_read) > max_linear
        # else:
        #     bad_data = self.image_stack > max_linear
        #
        # # normalize each image with the N-th frame to take out a linear slope
        # # normalized = self.image_stack / (self.image_stack[ref_frame_id] / ref_frame_id)
        #
        # # mask out all pixels above a saturation level
        # # normalized[bad_data] = numpy.NaN
        #
        # masked = self.image_stack.copy()
        # masked[bad_data] = numpy.NaN
        #
        # linearized_intensity = numpy.arange(masked.shape[0]).reshape((-1,1,1)) * (self.image_stack[ref_frame_id:ref_frame_id+1] / ref_frame_id)
        # print(linearized_intensity.shape)
        # pyfits.PrimaryHDU(data=linearized_intensity).writeto("linearized.fits", overwrite=True)
        #

        # now fit a each pixel

        parallel = False

        if (parallel):
            _iy,_ix = numpy.indices((masked.shape[1], masked.shape[2]))
            _ixy = numpy.dstack([_ix,_iy]).reshape((-1,2))
            print("ixy shape:", _ixy.shape)
            _ixy = _ixy[:25]
            # print(_ixy[:25])


            pool = multiprocessing.Pool(processes=2) #multiprocessing.cpu_count())

            __ixy = list(zip(_ixy[:,0], _ixy[:,1]))
            _masked = [masked[iy,ix] for (ix,iy) in __ixy]
            _linint = [linearized_intensity[iy,ix] for (ix,iy) in __ixy]
            # print(masked)
            # print(list(zip(_ixy[:,0], _ixy[:,1])))

            print("masked")
            print(_masked)
            print(len(_masked))
            print(len(list(_masked)))

            print("linint")
            print(_linint)
            print(len(list(_linint)))

            print("it")
            it = zip(_masked, _linint)
            print(it)
            print(len(list(it)))

            results_parallel = pool.starmap(
                self._fit_nonlinearity_pixel,
                iterable=zip(_masked, _linint),
            )
            pool.close()
            pool.join()

            # [masked[:, iy, ix], linearized_intensity[:, iy, ix] for [ix,iy] in _ixy],
            #     # iterable=zip(itertools.repeat(masked),
            #     #              itertools.repeat(linearized_intensity),
            #     #              _ixy[:, 0], _ixy[:, 1]),

            print("results_parallel=\n", results_parallel)
        else:
            cube_shape = self.image_stack.shape
            nonlinearity_fits_3d = numpy.zeros((3, cube_shape[1], cube_shape[2]))
            nonlinearity_fits_inverse = numpy.zeros((3, cube_shape[1], cube_shape[2]))
            for (ix, iy) in itertools.product(range(cube_shape[1]),
                                              range(cube_shape[2])):
            # for (ix, iy) in itertools.product(range(250,cube_shape[1],5),
            #                                   range(250,cube_shape[2],5)):
                if (iy == 0):
                    sys.stdout.write("\rWorking on column % 4d" % (ix))
                    sys.stdout.flush()

                    # print(ix, iy)
                # if ((ix % 100) == 0):
                #

                # sys.stdout.write("ix=% 4d  iy=% 4d\r" % (ix,iy))
                # sys.stdout.flush()
                # print(ix, iy)

                # extract data for this pixel
                raw_series = self.image_stack[:, iy, ix]
                series = raw_series - raw_series[0]

                diffs = numpy.pad(numpy.diff(raw_series), (1, 0))
                # print(diffs.shape, series.shape)

                # flag bad/saturated pixels
                max_diff = numpy.nanpercentile(diffs, 90)
                bad = (raw_series > 63000) | (diffs < 0.3 * max_diff)

                n_good = numpy.sum(~bad)
                if (n_good < 5):
                    continue

                avg_rate = series[10] / 10.
                # print(avg_rate)

                # perform initial fit to obtain a best-possible linear countrate
                integrations_count = numpy.arange(series.shape[0])
                pfit2 = self._fit_nonlinearity_pixel(integrations_count[~bad], series[~bad])
                best_fit_direct = pfit2[0] * integrations_count + \
                                  pfit2[1] * numpy.power(integrations_count, 2) + \
                                  pfit2[2] * numpy.power(integrations_count, 3)

                integrations_count = numpy.arange(series.shape[0])
                computed_countrate = integrations_count * pfit2[0]

                last_good_sample = numpy.max(integrations_count[~bad])

                # fit the non-linearity correction
                # print("Fitting pixel")
                pfit = self._fit_nonlinearity_pixel(series[~bad], computed_countrate[~bad])
                linearized = pfit[0] * numpy.power(series, 1) + \
                             pfit[1] * numpy.power(series, 2) + \
                             pfit[2] * numpy.power(series, 3)

                # _x = masked[:, iy, ix]
                # _y = linearized_intensity[:, iy, ix]
                # pfit = self._fit_nonlinearity_pixel(_x, _y)
                if (pfit is not None):
                    # nonlinearity_fits_3d[:, iy:iy+4, ix:ix+4] = pfit.reshape((-1,1,1))
                    nonlinearity_fits_3d[:, iy, ix] = pfit #.reshape((-1,1,1))

                pfit_inverse = self._fit_nonlinearity_pixel(computed_countrate[~bad], series[~bad])
                if (pfit is not None):
                    # nonlinearity_fits_inverse[:, iy:iy+4, ix:ix+4] = pfit_inverse.reshape((-1,1,1))
                    nonlinearity_fits_inverse[:, iy, ix] = pfit_inverse #.reshape((-1,1,1))

                if (make_plot):
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    # ax.scatter(_x, _y)
                    # ax.plot(_x, self._nonlinearity_fit_fct(pfit, _x))

                    ax.scatter(integrations_count, raw_series, s=4, label="raw data")
                    ax.scatter(integrations_count, series, s=8, label="zero subtracted")
                    ax.scatter(integrations_count[bad], series[bad], c='grey', marker='x', s=16, label='bad',
                               linewidth=1)
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
                    fig.tight_layout()
                    fig.suptitle("y = %+.4g*x %+.4g*x^2 %+.4g*x^3" % (pfit[0], pfit[1], pfit[2]))
                    fig.savefig("nonlin_plots/nonlin__x%04d__y%04d.png" % (ix, iy))
                    plt.close(fig)

            #break

        # return

        pyfits.PrimaryHDU(data=nonlinearity_fits_3d).writeto("nonlin3d.fits", overwrite=True)
        pyfits.PrimaryHDU(data=nonlinearity_fits_inverse).writeto("nonlin_inverse.fits", overwrite=True)
        return

    def read_nonlinearity_corrections(self, nonlin_fn):

        if (os.path.isfile(nonlin_fn)):
            try:
                hdulist = pyfits.open(nonlin_fn)
                print("Reading nonlinearity corrections from %s" % (nonlin_fn))
                nonlinearity_cube = hdulist[0].data
                print("CORR shape:", nonlinearity_cube.shape)
            except:
                return False

        self.nonlin_fn = nonlin_fn
        self.nonlinearity_cube = nonlinearity_cube

    def apply_nonlinearity_corrections(self, img_cube=None):

        if (self.nonlinearity_cube is None):
            print("No nonlinearity corrections loaded, skipping")
            return None

        # self.subtract_first_read()
        if (img_cube is None):
            img_cube = self.image_stack

        print("NONLIN: data=%s   corr=%s" % (str(img_cube.shape), str(self.nonlinearity_cube.shape)))

        linearized_cube = \
            self.nonlinearity_cube[0:1, :, :] * numpy.power(img_cube, 1) + \
            self.nonlinearity_cube[1:2, :, :] * numpy.power(img_cube, 2) + \
            self.nonlinearity_cube[2:3, :, :] * numpy.power(img_cube, 3)

        return linearized_cube

    def write_results(self, fn=None):
        if (fn is None):
            fn = os.path.join(self.basedir, self.filebase) + ".reduced.fits"
        hdulist = pyfits.HDUList([
            pyfits.PrimaryHDU(), #header=self.primary_header)
            pyfits.ImageHDU(data=self.weighted_mean, name="SCI"),
            pyfits.ImageHDU(data=self.noise_image, name='NOISE')
        ])
        print("Writing reduced results to %s" % (fn))
        hdulist.writeto(fn, overwrite=True)
        return

    def plot_pixel_curve(self, x, y, filebase=None,
                         cumulative=True, differential=False,
                         diff_vs_cum=False,
                         show_fits=False, show_errors=False):

        # self.subtract_first_read()
        counts = self.image_stack[:, y-1, x-1]
        zerolevel = self.image_stack[0, y-1,x-1]
        frame_number = numpy.arange(counts.shape[0])
        if (hasattr(self, "linearized_cube")):
            linearized_counts = self.linearized_cube[:, y - 1, x - 1]
        else:
            linearized_counts = counts
        phot_error = numpy.sqrt(counts - zerolevel)

        if (cumulative):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(frame_number, counts, s=4, label='raw')
            ax.scatter(frame_number, counts - zerolevel, s=4, label='raw, zerosub')
            ax.scatter(frame_number, linearized_counts, s=4, label='linearized, zerosub')
            ax.scatter(frame_number, linearized_counts + zerolevel, s=4, label='linearized')

            if (show_errors):
                ax.errorbar(frame_number, linearized_counts, yerr=phot_error)

            ax.axhline(y=63000, linestyle=":", color='grey')

            plot_fn = "cumpixelcurve_x%04d_y%04d.png" % (x, y)
            if (filebase is not None):
                plot_fn = filebase + plot_fn
            fig.savefig(plot_fn)
            plt.close(fig)

        if (differential):
            fig = plt.figure()
            ax = fig.add_subplot(111)

            diff = numpy.pad(numpy.diff(counts), (1,0), mode='constant', constant_values=0)

            ax.scatter(frame_number, diff, s=4, label='raw')
            #ax.scatter(frame_number, counts - zerolevel, s=4, label='raw, zerosub')
            #ax.scatter(frame_number, linearized_counts, s=4, label='linearized, zerosub')
            #ax.scatter(frame_number, linearized_counts + zerolevel, s=4, label='linearized')
            if (show_errors):
                ax.errorbar(frame_number, diff, yerr=phot_error)

            ax.axhline(y=0, linestyle=":", color='grey')

            plot_fn = "diffpixelcurve_x%04d_y%04d.png" % (x, y)
            if (filebase is not None):
                plot_fn = filebase + plot_fn
            fig.savefig(plot_fn)
            plt.close(fig)

        if (diff_vs_cum):
            fig = plt.figure()
            ax = fig.add_subplot(111)

            diff = numpy.pad(numpy.diff(counts), (1,0), mode='constant', constant_values=0)

            ax.scatter(counts, diff, s=4, label='raw')
            #ax.scatter(frame_number, counts - zerolevel, s=4, label='raw, zerosub')
            #ax.scatter(frame_number, linearized_counts, s=4, label='linearized, zerosub')
            #ax.scatter(frame_number, linearized_counts + zerolevel, s=4, label='linearized')

            if (show_errors):
                ax.errorbar(counts, diff, xerr=phot_error, yerr=phot_error)

            ax.axhline(y=0, linestyle=":", color='grey')
            ax.axvline(x=63000, linestyle=":", color='grey')
            plot_fn = "_diff_vs_cum__pixelcurve_x%04d_y%04d.png" % (x, y)
            if (filebase is not None):
                plot_fn = filebase + plot_fn
            fig.savefig(plot_fn)
            plt.close(fig)


    def dump_pixeldata(self, x, y, filebase=None, extras=None):

        print("dumping pixeldata for pixel @ %d / %d" % (x,y))

        _x = x-1
        _y = y-1
        frame_number = numpy.arange(self.image_stack.shape[0])
        raw_series = self.image_stack[:, _y, _x]
        linearized = self.linearized_cube[:, _y, _x]

        fn = "pixeldump_x%04d_y%04d.complete" % (x,y)
        if (filebase is not None):
            fn = filebase + fn

        extra_pixels = []
        if (extras is not None):
            try:
                for ex in extras:
                    extra_pixels.append(ex[:, _y, _x])
            except:
                pass

        numpy.savetxt(fn, numpy.array([
            frame_number, raw_series, linearized
            ]+extra_pixels).T
        )


if __name__ == "__main__":

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
    # cmdline.add_argument("healpix32", nargs="+", type=int,
    #                      help="list of input filenames")
    # cmdline.add_argument("--rerun", type=int, default=6,
    #                      help="rerun")
    cmdline.add_argument("--dumps", dest="write_dumps", default=False, action='store_true',
                         help="write intermediate process data [default: NO]")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    for fn in args.files:
        # fn = sys.argv[1]

        rss = RSS(fn, max_number_files=args.max_number_files)
        if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
            rss.read_nonlinearity_corrections(args.nonlinearity_fn)
        rss.reduce(write_dumps=args.write_dumps,
                   dark_fn=args.dark_fn)
        rss.write_results(fn="%s.%s.fits" % (rss.filebase, args.output_postfix))

        rss.plot_pixel_curve(818,1033)
        rss.plot_pixel_curve(1700,555)
        rss.plot_pixel_curve(250,1660)

        print("all done!")