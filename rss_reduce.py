#!/usr/bin/env python3

import sys
print(sys.path)

import logging
import multiprocessing
import os
import queue
import threading
import time
import glob
import datetime

import multiparallel_logging as mplog

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

import rss_filepicker

import astropy
print(astropy.__path__)


import rss_refpixel_calibrate

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

#
# Helper function for signal fitting
#

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


n_persistency_values = 8
def persistency_fit_pixel(differential_cube, linearized_cube, read_times, x, y):

    rate_series = differential_cube[:, y, x]
    linear_series = linearized_cube[:, y, x]

    # TODO: implement better noise model, accounting for read-noise and gain
    uncertainties = numpy.sqrt(linearized_cube[:, y, x])

    good4fit = numpy.isfinite(read_times) & \
               numpy.isfinite(rate_series) & \
               numpy.isfinite(uncertainties) & \
               (linear_series < 62000)

    read_time = read_times[good4fit]
    rate = rate_series[good4fit]
    uncert = uncertainties[good4fit]

    avg_rate = numpy.mean(rate)

    fallback_solution = [avg_rate, 0, 0]
    fallback_uncertainty = [0, 0, -1.]

    if (numpy.sum(good4fit) < 3):
        # if there's no good data we can't do any fitting
        return None,None,good4fit  # numpy.array(fallback_solution), numpy.array(fallback_uncertainty)  # assume perfect linearity

    # variables are: linear_rate, persistency_amplitude, persistency_timescale
    pinit = [numpy.min(rate), 2 * numpy.max(rate), 3.5]
    fit = scipy.optimize.leastsq(
        func=_persistency_plus_signal_fit_err_fct, x0=pinit,
        args=(read_time, rate, uncert),
        full_output=1
    )
    # print(fit)
    bestfit = fit[0]

    # Compute uncertainty on the shift and rotation
    if (fit[1] is not None):
        fit_uncert = numpy.sqrt(numpy.diag(fit[1]))
    else:
        fit_uncert = numpy.array([-99, -99., -99.])  # print(fit[1])

    special = False #(x>1025 & x<1050 & y>990 & y<1120)
    if (write_test_plot or special):
        fig = plt.figure()
        fig.suptitle("x=%d    y=%d" % (x,y))

        ax = fig.add_subplot(111)
        ax.scatter(read_times[good4fit], rate_series[good4fit], marker='o')
        ax.scatter(read_times[~good4fit], rate_series[~good4fit], marker='o', facecolors='none')
        ax.set_xlabel("Integration time")
        ax.set_ylabel("differential count")
        ax.set_yscale('log')
        # ax.set_xscale('log')
        #ax.set_ylim((numpy.max([1, numpy.min(rate_series[good4fit])]), 1.8 * numpy.max(rate_series[good4fit])))
        ystart = numpy.min([250, numpy.max([1, numpy.min(rate_series[good4fit])])])
        # ystart = numpy.max([10, numpy.min(rate_series[good4fit])])
        ax.set_ylim((ystart, 1.3 * numpy.max(rate_series[good4fit])))
        ax.plot(read_times, _persistency_plus_signal_fit_fct(bestfit, read_times))

        ax2 = ax.twinx()
        ax2.set_ylabel("linearized read counts")
        ax2.spines['right'].set_color('red')
        # ax2.spines['left'].set_color('blue')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        ax2.scatter(read_times[good4fit], linear_series[good4fit], c='red')
        ax2.scatter(read_times[~good4fit], linear_series[~good4fit], c='red', facecolors='none')
        ax2.axhline(y=62000, linestyle=":", color='red')

        ax.set_title("S(t) = %.1f + %.1f * exp(-t/%.3f)" % (bestfit[0], bestfit[1], bestfit[2]))

        plot_fn = "__debug_y=%04d_x=%04d.png" % (y,x)
        fig.savefig(plot_fn, bbox_inches='tight')
        plt.close(fig)

    return bestfit, fit_uncert, good4fit


def persistency_process_worker(
        row_queue,
        shmem_differential_cube, shmem_linearized_cube, shmem_persistency_fit,
        read_times,
        n_frames, nx=2048, ny=2048,
        name="Worker",
        write_test_plots=False,
    ):

    # make the shared memory available as numpy arrays
    linearized_cube = numpy.ndarray(
        shape=(n_frames,ny,nx), dtype=numpy.float32,
        buffer=shmem_linearized_cube.buf
    )
    differential_cube = numpy.ndarray(
        shape=(n_frames, ny, nx), dtype=numpy.float32,
        buffer=shmem_differential_cube.buf
    )
    persistency_fit = numpy.ndarray(
        shape=(n_persistency_values, ny, nx), dtype=numpy.float32,
        buffer=shmem_persistency_fit.buf,
    )

    logger = logging.getLogger("Persistency_%s" % (name))

    while (True):
        cmd = row_queue.get()
        if (cmd is None):
            break
        (row, full_fit_mask) = cmd

        logger.info("row % 4d: % 4d full fits (%d)" % (row, numpy.sum(full_fit_mask), nx))
        # print("row % 4d: % 4d full fits (%d)" % (row, numpy.sum(full_fit_mask), nx))

        linebuffer = numpy.full((n_persistency_values, nx), fill_value=numpy.NaN)

        for x in range(nx):
            # results = rss.fit_signal_with_persistency_singlepixel(
            #     x=x, y=row, debug=False, plot=False
            # )

            if (full_fit_mask[x]):
                # do a full fit for this pixel
                results = persistency_fit_pixel(
                    differential_cube=differential_cube,
                    linearized_cube=linearized_cube,
                    read_times=read_times,
                    x=x, y=row,
                )
                best_fit, fit_uncertainties, good4fit = results
                if (best_fit is not None):
                    linebuffer[0:3, x] = best_fit
                    linebuffer[3:6, x] = fit_uncertainties

                    integrated_persistency = \
                        best_fit[1] * best_fit[2] * (
                        numpy.exp(-read_times[1]/best_fit[2]) - numpy.exp(-numpy.nanmax(read_times)/best_fit[2]))
                    linebuffer[6, x] = integrated_persistency

                linebuffer[7, x] = numpy.sum(good4fit)

            else:
                # no need for a full fit, just calculate a simple slope
                diff_reads = differential_cube[:, row, x]
                raw_reads = linearized_cube[:, row, x]
                good_data = (raw_reads > 0) & (raw_reads < 62000)
                #print(diff_reads.shape)

                _median, _sigma = numpy.NaN, numpy.NaN
                for iter in range(3):
                    if (numpy.sum(good_data) < 1):
                        # no more good data left, so stick with what we had before
                        break
                    _stats = numpy.nanpercentile(diff_reads[good_data], [16,50,84])
                    _median = _stats[1]
                    _sigma = 0.5 * (_stats[2] - _stats[0])
                    outlier = (diff_reads > (_median + 3*_sigma)) | (diff_reads < (_median - 3*_sigma))

                    good_data[outlier] = False

                linebuffer[0,x] = _median
                linebuffer[3,x] = _sigma
                linebuffer[7,x] = numpy.sum(good_data)

                # linebuffer = numpy.array([_median, numpy.NaN, numpy.NaN,    # slope, persistency amp & time
                #                           _sigma, numpy.NaN, numpy.NaN,     # errors/uncertainties
                #                           numpy.NaN]                        # integrated persistency signal
                #                          )

        persistency_fit[:, row, :] = linebuffer

        # time.sleep(numpy.random.random(1) * 0.1)


    return


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
                 saturation_level=62000,
                 saturation_fraction=0.25, saturation_percentile=95,
                 use_reference_pixels=True,
                 mask_saturated_pixels=False):

        self.fn = fn
        self.filelist = []

        self.logger = logging.getLogger("RSS")

        self.use_reference_pixels = use_reference_pixels
        self.image_stack_initialized = False
        self.first_read_subtracted = False
        self.first_read = None
        self.first_header = None

        self.nonlin_fn = None
        self.nonlinearity_cube = None

        self.alloc_persistency = False

        # store values we may/will need during reduction
        self.max_number_files = -1 if max_number_files is None else max_number_files
        self.saturation_level = saturation_level
        self.saturation_fraction = saturation_fraction
        self.saturation_percentile = saturation_percentile
        self.mask_saturated_pixels = mask_saturated_pixels

        self.read_exposure_setup()

        self.get_full_filelist()

    def read_exposure_setup(self):
        if (self.fn is None):
            self.logger.critical("Unable to get exposure setup without valid input filename")

        # read the input file as reference file
        self.ref_hdulist = pyfits.open(self.fn)
        self.ref_header = self.ref_hdulist[0].header

        # image dimensions
        self.nx = self.ref_header['XSTOP'] - self.ref_header['XSTART'] + 1
        self.ny = self.ref_header['YSTOP'] - self.ref_header['YSTART'] + 1

        # readout settings
        self.n_groups = self.ref_header['NGROUPS']
        self.n_ramps = self.ref_header['NRAMPS']
        self.n_reads = self.ref_header['NREADS']
        self.n_outputs = self.ref_header['NOUTPUTS']
        self.gain = self.ref_header['GAIN']

        # exposure and other times
        self.exptime = self.ref_header['USEREXP'] / 1000.
        self.diff_exptime = self.exptime / self.n_groups


    def get_full_filelist(self):
        # get basedir
        fullpath = os.path.abspath(self.fn)
        print(fullpath)
        self.basedir, _fn = os.path.split(fullpath)
        self.filebase = ".".join(_fn.split(".")[:-2])
        # print(self.basedir, filebase)

        for _read in range(1,1000):
            filename = "%s.%d.fits" % (self.filebase, _read)
            full_filename = os.path.join(self.basedir, filename)
            if (os.path.isfile(full_filename)):
                self.filelist.append(full_filename)
            else:
                break

            # print(full_filename, os.path.isfile(full_filename))

        self.logger.debug("Loading filelist:\n"+"\n -- ".join(self.filelist))

        return

    def add_file(self, filename):
        return

    # def store_header_info(self, hdr):
    #     self.first_header = hdr
    #
    #     self.exptime = hdr['USEREXP'] / 1000.
    #     self.n_groups = hdr['NGROUPS']
    #     self.diff_exptime = self.exptime / self.n_groups

    def load_all_files(self, max_number_files=None, mask_saturated_pixels=False):

        if (max_number_files is None):
            max_number_files = self.max_number_files
        if (self.image_stack_initialized):
            self.logger.debug("stack already initialized, skipping repeat try")
            return

        self._image_stack = []

        # open all frames
        _filelist = self.filelist
        # if (max_number_files > 0):
        #     print("Limiting filelist to %d files" % (max_number_files))
        #     _filelist = _filelist[:max_number_files]

        # setup the data-cube to hold all the data
        if (max_number_files is not None and max_number_files > 0 and self.n_groups > max_number_files):
            self.n_groups = max_number_files
            print("Limiting input data to %d read-groups" % (max_number_files))

        self.image_stack_raw = numpy.full(
            (self.n_reads, self.n_groups, self.ny, self.nx),
            fill_value=numpy.NaN, dtype=numpy.float32)
        self.raw_read_times = numpy.full((self.n_reads, self.n_groups), fill_value=numpy.NaN)

        self.logger.debug("raw image cube dimensions: %s" % (str(self.image_stack_raw.shape)))

        # TODO: Add proper handling for combined Fowler and up-the-ramp sampling
        for fn in _filelist:
            hdulist = pyfits.open(fn)
            hdr = hdulist[0].header
            imgdata = hdulist[0].data.astype(numpy.float32)
            # hdulist.info()

            # mask all saturated pixels if requested
            if (mask_saturated_pixels):
                saturation_mask = (imgdata > self.saturation_level)
                print("masking out %d saturated pixels (%.1f)" % (
                    numpy.sum(saturation_mask), self.saturation_level))
                imgdata[saturation_mask] = numpy.Inf

            img_group = hdr['GROUP']
            img_read = hdr['READ']
            img_exptime = hdr['ACTEXP'] / 1000000. # convert time from raw microseconds to seconds
            self.logger.debug("FN=%s // grp=%d rd=%d exptime=%.4f" % (fn, img_group, img_read, img_exptime))

            if (max_number_files > 0 and img_group >= max_number_files):
                self.logger.debug("img-group > max-number-file --> skipping this file")
                continue

            self.raw_read_times[img_read-1, img_group-1] = img_exptime
            self.logger.debug("raw read times: %s" % (str(self.raw_read_times)))

            self.image_stack_raw[img_read-1, img_group-1, :, :] = imgdata



            # self._image_stack.append(imgdata)
            # if  (self.first_header is None):
            #     self.store_header_info(hdulist[0].header)

            # break

        # calculate the initial image stack
        self.logger.info("#groups=%d // #ramps=%d // #reads=%d" % (self.n_groups, self.n_ramps, self.n_reads))
        if (self.n_groups == 1 and self.n_reads >= 1):
            # this is a plain fowler mode, so calculate pair-wise differences
            self.logger.info("Using fowler-sampling strategy")
            self.image_stack = numpy.diff(self.image_stack_raw, axis=0)
            print("@@@@@@@@@@@@@", self.image_stack.shape)
            self.read_times = numpy.nanmean(self.raw_read_times, axis=0)

        elif (self.n_groups > 1):
            # up-the-ramp sampling.
            self.logger.info("Using up-the-ramp strategy")
            self.image_stack = numpy.nanmean(self.image_stack_raw, axis=0)
            print(self.raw_read_times)
            self.read_times = numpy.nanmean(self.raw_read_times, axis=0)

        else:
            self.logger.critical("No idea what's going here and what to do with this data - HELP!!!!")
            return


#         self.image_stack = numpy.array(self._image_stack, dtype=numpy.float32)

        # self.image_stack = numpy.diff(self.image_stack_raw, axis=0)

        self.logger.debug("stack before/after: %s --> %s" % (str(self.image_stack_raw.shape), str(self.image_stack.shape)))

        self.logger.info("read-times: %s" % (str(self.read_times)))

        self.logger.debug("stack shape: %s" % (str(self.image_stack.shape)))

        # delete raw stack to clean up memory
        del self.image_stack_raw

        self.image_stack_initialized = True

    def reduce(self, dark_fn=None, write_dumps=False, mask_bad_data=None, mask_saturated_pixels=False):

        self.load_all_files(mask_saturated_pixels=mask_saturated_pixels)

        if (self.use_reference_pixels):
            self.logger.info("Applying reference pixel corrections")
            reset_frame_subtracted = self.image_stack.copy()
            for frame_id in range(self.image_stack.shape[0]):
                reference_pixel_correction = rss_refpixel_calibrate.reference_pixels_to_background_correction(
                    self.image_stack[frame_id]
                )
                reset_frame_subtracted[frame_id] -= reference_pixel_correction
        else:
            self.logger.info("Subtracting first read from stack")
            # self.subtract_first_read()
            # apply first-read subtraction
            self.reset_frame = self.image_stack[0]
            reset_frame_subtracted = self.image_stack - self.reset_frame

        # apply any necessary corrections for nonlinearity and other things
        self.logger.info("Applying non-linearity corrections")
        linearized = self.apply_nonlinearity_corrections(reset_frame_subtracted)
        # print("linearized = ", linearized)
        if (linearized is None):
            self.logger.warning("No linearized data found, using raw data instead")
            linearized = reset_frame_subtracted

        # prepare shared memory to receive the linearized data cube
        self.logger.debug("Allocating shared memory for linearized cube")
        self.shmem_linearized_cube = multiprocessing.shared_memory.SharedMemory(
            create=True, size=linearized.nbytes
        )
        self.linearized_cube = numpy.ndarray(
            shape=linearized.shape, dtype=numpy.float32,
            buffer=self.shmem_linearized_cube.buf
        )
        self.linearized_cube[:,:,:] = linearized[:,:,:]
        self.logger.debug("linearized cube initialized")

        dark_cube = numpy.zeros_like(linearized)
        if (dark_fn is None):
            self.logger.warning("No dark correction requested, skipping")
        elif (not os.path.isfile(dark_fn)):
            self.logger.warning("Dark requested (%s) but not found" % (dark_fn))
        else:
            try:
                # load dark-rate image
                self.logger.info("Loading dark-corrections from %s" % (dark_fn))
                dark_hdu = pyfits.open(dark_fn)
                dark = dark_hdu['DARKRATE'].data

                # perform a dark subtraction;
                # dark-current = rate [in cts/sec] * frame-# * exposure-time per frame [in sec]
                dark_cube = (numpy.arange(linearized.shape[0], dtype=numpy.float).reshape((-1,1,1)) + 1) \
                            * self.diff_exptime \
                            * dark.reshape((1, dark.shape[0], dark.shape[1]))
                self.logger.debug("shape of dark cube: %s" % (dark_cube.shape))
                self.linearized_cube -= dark_cube
            except Exception as e:
                self.logger.error("error during dark subtraction:\n%s" % (str(e)))

        # allocate shared memory for the differential stack and calculate from
        # the linearized cube
        # calculate differential stack
        self.logger.debug("Allocating shared memory for differential cube")
        self.shmem_differential_cube = multiprocessing.shared_memory.SharedMemory(
            create=True, size=self.linearized_cube.nbytes
        )
        self.differential_cube = numpy.ndarray(
            shape=self.linearized_cube.shape, dtype=numpy.float32,
            buffer=self.shmem_differential_cube.buf,
        )
        self.logger.debug("differential cube allocated")

        self.logger.debug("calculating differential cube")
        self.differential_cube[:, :, :] = numpy.pad(
            numpy.diff(linearized, axis=0), ((1,0),(0,0),(0,0))
        ) / self.read_times.reshape((-1,1,1))
        self.logger.debug("diff stack: %s" % (str(self.differential_cube.shape)))

        # mask out all saturated and/or otherwise bad samples
        max_count_rates = -1000 # TODO: FIX THIS numpy.nanpercentile(self.differential_stack, q=self.saturation_percentile, axis=0)
        # print("max counrates:", max_count_rates.shape)

        # TODO: implement full iterative outlier rejection here
        if (mask_bad_data is None):
            mask_bad_data = self.mask_BAD_DARK | self.mask_SATURATED | self.mask_LOW_RATE | self.mask_NEGATIVE
        self.logger.info("Identifying bad/dead/saturated/negative pixels (0x%02x)" % (mask_bad_data))
        bad_data = numpy.zeros_like(self.image_stack, dtype=numpy.bool)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_SATURATED) > 0):
            bad_data = bad_data | (self.image_stack > self.saturation_level)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_LOW_RATE) > 0):
            bad_data = bad_data | (self.differential_cube < self.saturation_fraction * max_count_rates)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_BAD_DARK) > 0):
            bad_data = bad_data | (dark_cube >= linearized)
        if (mask_bad_data is not None and (mask_bad_data & self.mask_NEGATIVE) > 0):
            bad_data = bad_data | (linearized < 0)

        self.bad_data_mask = bad_data

        # bad_data = (self.image_stack > self.saturation_level) | \
        #            (self.differential_stack < self.saturation_fraction*max_count_rates) | \
        #            (dark_cube >= linearized) | \
        #            (linearized < 0)

        self.logger.info("Cleaning image cube")
        self.clean_stack = self.differential_cube.copy()
        self.clean_stack[bad_data] = numpy.NaN
        self.clean_stack[0, :, :] = numpy.NaN # mask out the first slice, which is just padding

        # calculate a average countrate image
        self.logger.info("calculating final image from stack")
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
            self.logger.info("Writing all dumps")
            bn = self.filebase + "__"
            self.logger.debug("Dump-file basename: %s" % (bn))
            pyfits.PrimaryHDU(data=self.image_stack).writeto(bn+"stack_raw.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.image_stack).writeto(bn+"stack_zerosub.fits", overwrite=True)
            pyfits.PrimaryHDU(data=linearized).writeto(bn+"stack_linearized.fits", overwrite=True)
            self.logger.debug("writing darkcube")
            pyfits.PrimaryHDU(data=dark_cube).writeto(bn+"stack_darkcube.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.differential_cube).writeto(bn + "stack_diff.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.clean_stack).writeto(bn+"stack_clean.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=ratios).writeto("stack_ratios.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=self.reduced_image_plain).writeto("final_image.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=image7).writeto("final_image7.fits", overwrite=True)
            # pyfits.PrimaryHDU(data=max_count_rates).writeto("max_count_rates.fits", overwrite=True)
            pyfits.PrimaryHDU(data=bad_data.astype(numpy.int)).writeto(bn+"bad_data.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.weighted_mean).writeto(bn+"final_image_weighted.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.inv_noise).writeto(bn+"final_inv_noise.fits", overwrite=True)

            n_good_pixels = numpy.sum(~bad_data, axis=0)
            print("#goodpixels", n_good_pixels.shape)
            pyfits.PrimaryHDU(data=n_good_pixels).writeto(bn+"n_good_pixels.fits", overwrite=True)

        return

    def subtract_first_read(self):
        if (self.first_read_subtracted):
            self.logger.debug("First read already subtracted, skipping")
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
                self.logger.info("Reading nonlinearity corrections from %s" % (nonlin_fn))
                nonlinearity_cube = hdulist[0].data
                self.logger.debug("CORR shape: %s" % (nonlinearity_cube.shape))
            except:
                return False

        self.nonlin_fn = nonlin_fn
        self.nonlinearity_cube = nonlinearity_cube

    def apply_nonlinearity_corrections(self, img_cube=None):

        if (self.nonlinearity_cube is None):
            self.logger.warning("No nonlinearity corrections loaded, skipping")
            return None

        # self.subtract_first_read()
        if (img_cube is None):
            img_cube = self.image_stack

        self.logger.debug("NONLIN: data=%s   corr=%s" % (str(img_cube.shape), str(self.nonlinearity_cube.shape)))

        linearized_cube = \
            self.nonlinearity_cube[0:1, :, :] * numpy.power(img_cube, 1) + \
            self.nonlinearity_cube[1:2, :, :] * numpy.power(img_cube, 2) + \
            self.nonlinearity_cube[2:3, :, :] * numpy.power(img_cube, 3)

        return linearized_cube


    def write_results(self, fn=None, flat4salt=False):
        if (fn is None):
            fn = os.path.join(self.basedir, self.filebase) + ".reduced.fits"

        # collect all output results
        _list = [pyfits.PrimaryHDU(header=self.ref_header)]
        if (flat4salt):
            # only write the reduced frame and nothing else
            try:
                img = self.persistency_fit_global[:, :, 0]
                hdu = pyfits.ImageHDU(data=img, name="SCI")
                hdr = hdu.header
                hdr['FIT_PERS'] = (True, "true persistency results")
            except:
                img = self.weighted_mean
                hdu = pyfits.ImageHDU(data=img, name="SCI")
                hdr = hdu.header
                hdr['FIT_PERS'] = (False, "true persistency results")
            _list.append(hdu)
        else:
            _list.extend([
                pyfits.ImageHDU(data=self.weighted_mean, name="SCI"),
                pyfits.ImageHDU(data=self.noise_image, name='NOISE')
            ])
            try:
                for i,extname in enumerate([
                    'PERS.SIGNAL', 'PERS.AMP', 'PERS.TAU',
                    'PERS.ERR.SIGNAL', 'PERS.ERR.AMP', 'PERS.ERR.TAU',
                    'PERS.INTEGRATED', 'PERS.NSAMPLES']):
                    _list.append(
                        pyfits.ImageHDU(
                            data=rss.persistency_fit_global[i, :, :],
                            name=extname)
                    )
            except:
                pass

        hdulist = pyfits.HDUList(_list)
        self.logger.info("Writing reduced results to %s" % (fn))
        hdulist.writeto(fn, overwrite=True)
        return

    def _alloc_persistency(self):
        # allocate a datacube for the persistency fit results in shared memory
        # to make it read- and write-accessible from across all worker processes
        self.shmem_persistency_fit_global = multiprocessing.shared_memory.SharedMemory(
            create=True, size=n_persistency_values*self.nx*self.ny*4,
        )
        self.persistency_fit_global = numpy.ndarray(
            shape=(n_persistency_values, self.ny, self.ny), dtype=numpy.float32,
            buffer=self.shmem_persistency_fit_global.buf,
        )
        self.persistency_fit_global[:,:,:] = numpy.NaN
        self.alloc_persistency = True

    def fit_signal_with_persistency(
            self,
            n_workers=0,
            previous_frame=None,
            write_test_plots=False
    ):

        # by default use all existing CPU cores for processing
        if (n_workers <= 0):
            n_workers = multiprocessing.cpu_count()

        if (not self.alloc_persistency):
            self._alloc_persistency()

        # prepare and fill job-queue - split work into chunks of individual lines
        row_queue = multiprocessing.JoinableQueue()

        print("xxx")
        if (previous_frame is not None and os.path.isfile(previous_frame)):
            self.logger.info("Using previous frame (%s) to speed up persistency correction" % (previous_frame))
            prev_hdu = pyfits.open(previous_frame)
            prev_img = prev_hdu[0].data
            threshold = 58000
            need_full_persistency_fit = prev_img > threshold
            for y in numpy.arange(self.ny):
                row_queue.put((y, need_full_persistency_fit[y, :]))
        else:
            for y in numpy.arange(self.ny):
                row_queue.put((y, numpy.ones((self.nx), dtype=bool)))

        # setup threads, fill queue with stop signals, but don't start threads just yet
        fit_threads = []
        for n in range(n_workers):
            t = multiprocessing.Process(
                target=persistency_process_worker,
                kwargs=dict(
                    row_queue=row_queue,
                    shmem_differential_cube=self.shmem_differential_cube,
                    shmem_linearized_cube=self.shmem_linearized_cube,
                    shmem_persistency_fit=self.shmem_persistency_fit_global,
                    read_times=self.read_times,
                    n_frames=self.linearized_cube.shape[0],
                    nx=self.nx, ny=self.ny,
                    name="FitWorker_%02d" % (n+1),
                    write_test_plots=write_test_plots,
                )
            )
            t.daemon = True
            fit_threads.append(t)
            row_queue.put(None)

        # once all threads are setup, start them
        for t in fit_threads:
            t.start()

        # and then wait until they are all done with their work
        for t in fit_threads:
            t.join()

        # for now write the fit output
        self.logger.debug("dumping fit results")
        out_tmp = pyfits.PrimaryHDU(data=self.persistency_fit_global)
        out_tmp.writeto("persistency_fit_dump.fits", overwrite=True)
        return

    def fit_signal_with_persistency_singlepixel(self, x, y, debug=False, plot=False):

        _x = x - 1
        _y = y - 1

        bad_data = self.bad_data_mask[:, _y, _x]
        rate_series = self.differential_cube[:, _y, _x]

        # TODO: implement better noise model, accounting for read-noise and gain
        uncertainties = numpy.sqrt(self.image_stack[:, _y, _x])

        if (debug):
            numpy.savetxt(
                "persistency_dump_%04dx%04d.txt" % (x, y),
                numpy.array([self.read_times, rate_series, uncertainties, bad_data]).T,
            )

        good4fit = numpy.isfinite(self.read_times) & \
                   numpy.isfinite(rate_series) & \
                   numpy.isfinite(uncertainties) & \
                   ~bad_data
        read_time = self.read_times[good4fit]
        rate = rate_series[good4fit]
        uncert = uncertainties[good4fit]

        avg_rate = numpy.mean(rate)

        fallback_solution = [avg_rate, 0, 0]
        fallback_uncertainty = [0,0,-1.]

        if (numpy.sum(good4fit) < 5):
            # if there's no good data we can't do any fitting
            return None #numpy.array(fallback_solution), numpy.array(fallback_uncertainty)  # assume perfect linearity

        # variables are: linear_rate, persistency_amplitude, persistency_timescale
        pinit = [numpy.min(rate), 2 * numpy.max(rate), 3.5]
        fit = scipy.optimize.leastsq(
            func=_persistency_plus_signal_fit_err_fct, x0=pinit,
            args=(read_time, rate, uncert),
            full_output=1
        )
        # print(fit)
        bestfit = fit[0]

        # Compute uncertainty on the shift and rotation
        if (fit[1] is not None):
            fit_uncert = numpy.sqrt(numpy.diag(fit[1]))
        else:
            fit_uncert = numpy.array([-99, -99., -99.]) #print(fit[1])

        #print(numpy.diag(fit[1]))


        # bestfit = fit_persistency_plus_signal_pixel(
        #     rss.read_times[~bad_data], rate_series[~bad_data], uncertainties[~bad_data]
        # )
        # print("BESTFIT:", x, y, bestfit, "   /// uncertainties: ", fit_uncert)

        if (plot):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(self.read_times[~bad_data], rate_series[~bad_data], s=8, c='blue')
            ax.scatter(self.read_times[bad_data], rate_series[bad_data], s=4, c='grey')
            timeaxis = numpy.linspace(0, numpy.nanmax(self.read_times), 250)
            # print("read-times = \n", timeaxis)
            modely = _persistency_plus_signal_fit_fct(bestfit, timeaxis)
            # print("best-fit:", bestfit)
            # print("model-fit = \n", modely)

            ax.plot(timeaxis, modely)
            plot_fn = "%s____persistency_plus_signal__%04d-%04d.png" % (self.filebase, x, y)
            fig.suptitle("F = %.0f + %.0f x exp(-t/%.3fs)" % (bestfit[0], bestfit[1], bestfit[2]))
            ax.set_xlabel("integration time [seconds]")
            ax.set_ylabel("flux above read #0 [counts]")
            fig.savefig(plot_fn, dpi=200)
            plt.close(fig)

        return bestfit, fit_uncert

    def plot_pixel_curve(self, x, y, filebase=None,
                         cumulative=True, differential=False,
                         diff_vs_cum=False,
                         show_fits=False, show_errors=False,
                         show_plot=False):

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
            ax.legend()

            fig.suptitle("%s :: x=%d y=%d" % (self.filebase,x,y))

            plot_fn = "_pixelcurve_x%04d_y%04d____cum_vs_read.png" % (x, y)
            if (filebase is not None):
                plot_fn = filebase + plot_fn
            fig.savefig(plot_fn)

            if (show_plot):
                fig.show()
            else:
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
            ax.legend()

            fig.suptitle("%s :: x=%d y=%d" % (self.filebase,x,y))

            plot_fn = "_pixelcurve_x%04d_y%04d____diff_vs_read.png" % (x, y)
            if (filebase is not None):
                plot_fn = filebase + plot_fn
            fig.savefig(plot_fn)

            if (show_plot):
                fig.show()
            else:
                plt.close(fig)
            # plt.close(fig)

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

            fig.suptitle("%s :: x=%d y=%d" % (self.filebase,x,y))
            ax.legend()

            ax.axhline(y=0, linestyle=":", color='grey')
            ax.axvline(x=63000, linestyle=":", color='grey')
            plot_fn = "_pixelcurve_x%04d_y%04d____diff_vs_cum.png" % (x, y)
            if (filebase is not None):
                plot_fn = filebase + plot_fn
            fig.savefig(plot_fn)

            if (show_plot):
                fig.show()
            else:
                plt.close(fig)
            # plt.close(fig)


    def dump_pixeldata(self, x, y, filebase=None, extras=None):

        self.logger.debug("dumping pixeldata for pixel @ %d / %d" % (x,y))

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


    def _parallel_worker(self, job_queue, result_queue, execute_function, is_in_class=True):
        self.logger.debug("Worker has started")
        while (True):
            job = job_queue.get()
            if (job is None):
                job_queue.task_done()
                self.logger.debug("Worker shutting down")
                break

            [x, y] = job
            if ((y % 10) == 0 and x == 0): print(job)

            if (is_in_class):
                result = execute_function(x, y)
            else:
                result = execute_function(self, x, y)

            result_queue.put((x, y, result))

            job_queue.task_done()


    def parallel_fitter(self, xrange=None, yrange=None,
                        execute_function=None, return_dim=1, is_in_class=True,
                        n_workers=None):

        if (xrange is None):
            x1 = 0
            x2 = self.image_stack.shape[2]
        else:
            [x1,x2] = xrange

        if (yrange is None):
            y1 = 0
            y2 = self.image_stack.shape[1]
        else:
            [y1,y2] = yrange

        # prepare return results
        return_results = numpy.full((return_dim, y2-y1, x2-x1), fill_value=numpy.NaN)

        # generate the coordinates where we need to fit/process data
        # iy, ix = numpy.indices((rss.image_stack.shape[1], rss.image_stack.shape[2]))
        iy, ix = numpy.indices((y2-y1, x2-x1))
        iy += y1
        ix += x1
        # self.logger.debug(iy.shape)

        # prepare work and results queue for data exchange with the workers
        job_queue = multiprocessing.JoinableQueue()
        result_queue = multiprocessing.Queue()
        ixy = numpy.dstack([ix, iy]).reshape((-1, 2))
        # self.logger.debug(ixy.shape)
        # self.logger.debug(ixy[:10])

        # prepare and start the workers
        self.logger.debug("Creating workers")
        worker_processes = []
        if (n_workers is None):
            n_workers = multiprocessing.cpu_count()
        for n in range(n_workers):
            p = multiprocessing.Process(
                target=self._parallel_worker,
                kwargs=dict(#self=self,
                            job_queue=job_queue,
                            result_queue=result_queue,
                            execute_function=execute_function,
                            is_in_class=is_in_class
                            )
            )
            p.daemon = True
            p.start()
            worker_processes.append(p)

        # prepare jobqueue
        self.logger.debug("preparing jobs")
        n_jobs = 0
        for _xy in ixy:
            job_queue.put(_xy)
            n_jobs += 1

        # add termination signals to the work-queue so workers know when to shut down
        for p in worker_processes:
            job_queue.put(None)

        # wait for completion
        self.logger.debug("Waiting for completion")
        job_queue.join()

        # receive all the results back from the workers
        for j in range(n_jobs):
            result = result_queue.get()
            (x,y,value) = result
            return_results[:, y-y1,x-x1] = value

        return return_results



    def fit_2component_persistency_plus_signal(self, x, y):

        _x = x-1
        _y = y-1

        raw_data = self.image_stack[:, _y, _x]
        bad_data = raw_data > self.saturation_level
        # bad_data = self.bad_data_mask[:, _y, _x]
        if (hasattr(self, 'linearized_cube')):
            series = self.linearized_cube[:, _y, _x]
        else:
            series = self.image_stack[:,_y,_x]

        if (numpy.sum(~bad_data) < 5):
            return [numpy.NaN, numpy.NaN, numpy.NaN]

        pinit = [1., 0.] #, 1.]

        readout_times = numpy.arange(series.shape[0], dtype=numpy.float) * self.diff_exptime
        img_time = readout_times[~bad_data]
        img_flux = series[~bad_data]

        fit = scipy.optimize.leastsq(
            func=_persistency_plus_signal_fit_err_fct, x0=pinit,
            args=(img_time, img_flux),
            full_output=1
        )
        # print(fit)
        pfit = fit[0]

        # bestfit = _fit_persistency_plus_signal_pixel(
        #     integ_exp_time[~bad_data], series[~bad_data])

        return pfit #self.image_stack[1, y-1, x-1]

    def load_precalculated_results(self, weighted_image_fn=None, persistency_fit_fn=None):

        if (weighted_image_fn is not None and os.path.isfile(weighted_image_fn)):
            self.logger.warning("Loading weighted results from file isn't implemented yet")
            pass

        if (persistency_fit_fn is not None and os.path.isfile(persistency_fit_fn)):
            self.logger.info("Loading canned persistency results from %s" % (persistency_fit_fn))
            # make sure we have compatible memory
            if (not self.alloc_persistency):
                self._alloc_persistency()

            # read FITS file and copy the image into the allocated memory buffer
            self.logger.info("Loading pre-calculated persistency results [%s]" % (
                persistency_fit_fn)
            )
            hdulist = pyfits.open(persistency_fit_fn)
            self.persistency_fit_global[:,:,:] = hdulist[0].data[:,:,:]
        else:
            self.logger.warning("Unable to load previous persistency results (%s)" % (persistency_fit_fn))

            pass

    def find_previous_exposure(self, search_dir):

        self.logger.debug("Finding previous exposure (%s)" % (search_dir))

        dir_index = rss_filepicker.PreviousFilePicker(search_dir)
        prior_fn, delta_seconds = dir_index.find_closest(self.ref_header)

        if (prior_fn is None):
            self.logger.warn("Unable to find prior frame as persistency reference")
        else:
            self.logger.info("Found %s, taken %.2f seconds prior to %s" % (prior_fn, delta_seconds, self.ref_header['FILE']))

        return prior_fn, delta_seconds


    def __del__(self):
        # self.logger.debug("Running destructor and cleaning up shared memory")
        # clean up shared memory
        try:
            self.shmem_linearized_cube.close()
            self.shmem_linearized_cube.unlink()
        except:
            pass

        try:
            self.shmem_differential_cube.close()
            self.shmem_differential_cube.unlink()
        except:
            pass

        try:
            self.shmem_persistency_fit_global.close()
            self.shmem_persistency_fit_global.unlink()
        except:
            pass


if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("RunAnalysis")

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

    # cmdline.add_argument("healpix32", nargs="+", type=int,
    #                      help="list of input filenames")
    # cmdline.add_argument("--rerun", type=int, default=6,
    #                      help="rerun")
    cmdline.add_argument("--dumps", dest="write_dumps", default=False, action='store_true',
                         help="write intermediate process data [default: NO]")
    cmdline.add_argument("--debugpngs", dest="write_debug_pngs", default=False, action='store_true',
                         help="generate debug plots for all pixels with persistency [default: NO]")
    cmdline.add_argument("--refpixel", dest="use_ref_pixels", default=False, action='store_true',
                         help="use reference pixels [default: NO]")
    cmdline.add_argument("--flat4salt", dest="write_flat_for_salt", default=False, action='store_true',
                         help="write a flat, 1-extension FITS file for SALT")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    for fn in args.files:
        # fn = sys.argv[1]

        rss = RSS(fn, max_number_files=args.max_number_files,
                   use_reference_pixels=args.use_ref_pixels)
        if (args.nonlinearity_fn is not None and os.path.isfile(args.nonlinearity_fn)):
            rss.read_nonlinearity_corrections(args.nonlinearity_fn)
        rss.reduce(write_dumps=args.write_dumps,
                   dark_fn=args.dark_fn)

        persistency_options = args.persistency_mode.split(":")
        persistency_mode = persistency_options[0]
        have_persistency_results = False
        if (persistency_mode == "none"):
            logger.info("Nothing to do")
        elif (persistency_mode == "full"):
            logger.info("Calculating persistency results for all pixels")
            rss.fit_signal_with_persistency(previous_frame=None)
            have_persistency_results = True
        elif (persistency_mode == "best"):
            logger.info("Using on-demand persistency fitting")
            if(len(persistency_options) < 2):
                logger.info("Insufficient information, defaulting to running on all pixels")
                rss.fit_signal_with_persistency(previous_frame=None)
                have_persistency_results = True
            else:
                opt = persistency_options[1]
                if (os.path.isfile(opt)):
                    logger.info("Using optimized persistency mode (ref-fn: %s)" % (opt))
                    rss.fit_signal_with_persistency(previous_frame=opt)
                    have_persistency_results = True
                elif (os.path.isdir(opt)):
                    logger.info("Searching for optimal reference frame in --> %s <--" % (opt))
                    xxx = rss.find_previous_exposure(opt)  #find_previous_frame(rss.ref_header, opt)
                    # print(xxx)
                    ref_fn, delta_t = xxx #rss.find_previous_exposure(opt)  #find_previous_frame(rss.ref_header, opt)
                    if (ref_fn is not None):
                        logger.info("Using optimized persistency mode using automatic ref-fn: %s (Dt=%.3f)" % (ref_fn, delta_t))
                        rss.fit_signal_with_persistency(previous_frame=ref_fn,
                                                        write_test_plots=args.write_debug_pngs)
                        have_persistency_results = True
                    else:
                        logger.warning("No previous frame found, skipping persistency fit")
                        # rss.fit_signal_with_persistency(previous_frame=ref_fn)
                        # have_persistency_results = True
                else:
                    logger.info("Unknown option to best mode (found: %s), skipping persistency modeling" % (opt))
        else:
            logger.info("Unknown persistency request (%s)" % (persistency_mode))

        if (have_persistency_results):
            out_tmp = pyfits.PrimaryHDU(data=rss.persistency_fit_global)
            fit_fn = "%s.%s.persistencyfit.fits" % (rss.filebase, args.output_postfix)
            logger.info("Writing persistency fit to %s ..." % (fit_fn))
            out_tmp.writeto(fit_fn, overwrite=True)

        red_fn = "%s.%s.fits" % (rss.filebase, args.output_postfix)
        logger.info("Writing reduction results to %s" % (red_fn))
        rss.write_results(fn=red_fn, flat4salt=args.write_flat_for_salt)

        # rss.plot_pixel_curve(818,1033)
        # rss.plot_pixel_curve(1700,555)
        # rss.plot_pixel_curve(505,1660)

        del rss
        print("all done!")