#!/usr/bin/env python3
import queue

import sys
# print(sys.path)

import logging
import os

import multiparlog as mplog

import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import itertools
import multiprocessing
import multiprocessing.shared_memory
import argparse
import time

from astropy import log
log.setLevel('ERROR')

import warnings
warnings.filterwarnings('ignore')

from .provenance import DataProvenance
from .refpixel_calibrate import  reference_pixels_to_background_correction

from .nirwals_urg_algorithms import *

import nirwals.data

import astropy


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

persistency_values = [
    'PERS.SIGNAL', 'PERS.AMP', 'PERS.TAU', 'PERS2.AMP', 'PERS2.TAU',
    'PERS.SIGNAL.ERR', 'PERS.AMP.ERR', 'PERS.TAU.ERR', 'PERS2.AMP.ERR', 'PERS2.TAU.ERR',
    'PERS.INTEGRATED', 'N_PIXELS', 'RAMP.OFFSET', 'SIGNAL_FRACTION'
]
n_persistency_values = len(persistency_values)
n_pars = 5

def persistency_fit_pixel(differential_cube, linearized_cube, read_times,
                          good_data,
                          x, y, write_test_plot=False):

    rate_series = differential_cube[:, y, x]
    linear_series = linearized_cube[:, y, x]

    # TODO: implement better noise model, accounting for read-noise and gain
    uncertainties = numpy.sqrt(linearized_cube[:, y, x])

    good4fit = numpy.isfinite(read_times) & \
               numpy.isfinite(rate_series) & \
               numpy.isfinite(uncertainties) & \
               (linear_series < 55000)
    if (good_data is not None):
        good4fit = good4fit & good_data

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

    # work out best initial guesses
    # pinit = [numpy.min(rate), 2 * numpy.max(rate), 3.5]
    # rate: minimum encountered, but minimum 0, max 60000
    rate_guess = numpy.max([0., numpy.min([60000, numpy.nanmin(rate)])])
    # persistency: max-rate minus best-guess for signal rate
    pers_guess = numpy.max([0, numpy.min([60e3, (numpy.nanmax(rate) - rate_guess)])])
    # tau: typical value found in data fits
    tau_guess = 1.4
    pinit = [rate_guess, pers_guess, tau_guess]

    # boundary conditions for exponential fit parameters:
    # signal:   -10 ... Inf counts/sec (small negative allowed for noise)
    # persisntency amplite: 0 .. 65K counts/sec
    # timescale tau: 0.2 .. 100 seconds
    pinit = numpy.array([rate_guess, pers_guess, 2.0, 0, 50.])
    try:
        fit_results = scipy.optimize.least_squares(
            fun=_persistency_plus_signal_fit_err_fct2,
            x0=pinit,
            bounds=([0., 0., 1., 0., 5.], [numpy.Inf, 65e3, 4, 65e3, 2000]),
            kwargs=dict(read_time=read_time, rate=rate,
                        uncert=uncert),
        )

        # fit_results = scipy.optimize.least_squares(
        #     fun=_persistency_plus_signal_fit_err_fct,
        #     x0=pinit,
        #     bounds=([-100, 0, 0.2], [numpy.Inf, 65e3, 1000.]),
        #     kwargs=dict(read_time=read_time, rate=rate, uncert=uncert),
        # )
        bestfit = fit_results.x
        bounds_limited_mask = fit_results.active_mask
        fit_successful = fit_results.success
    except ValueError as e:
        bestfit = [numpy.NaN,numpy.NaN,numpy.NaN]
        bounds_limited_mask = False
        fit_successful = False


    # fit = scipy.optimize.leastsq(
    #     func=_persistency_plus_signal_fit_err_fct, x0=pinit,
    #     args=(read_time, rate, uncert),
    #     full_output=1
    # )
    # # print(fit)
    # bestfit = fit[0]

    # Compute uncertainty on the shift and rotation
    if (fit_successful):
        # TODO: Figure out uncertainties here
        fit_uncert = numpy.zeros_like(pinit)
        #[0,0,0] #numpy.sqrt(numpy.diag(fit[1]))
    else:
        fit_uncert = numpy.full_like(pinit, fill_value=-99)
        #([-99, -99., -99.])  # print(fit[1])

    special = False #(x>1025 & x<1050 & y>990 & y<1120)
    # special = (x>1650 and x<1700 and y>400 and y<450)
    # special = (x >380 and x<500 and y > 200 and y < 250)
    if (write_test_plot or special):
        fig = plt.figure()
        fig.suptitle("x=%d    y=%d" % (x,y))

        ax = fig.add_subplot(111)
        ax.scatter(read_times[good4fit], rate_series[good4fit], marker='o', s=2)
        ax.scatter(read_times[~good4fit], rate_series[~good4fit], marker='o', facecolors='none', s=2, alpha=0.5)
        ax.set_xlabel("Integration time")
        ax.set_ylabel("differential count")
        # ax.set_yscale('log')
        # ax.set_xscale('log')
        #ax.set_ylim((numpy.max([1, numpy.min(rate_series[good4fit])]), 1.8 * numpy.max(rate_series[good4fit])))
        ystart = numpy.min([250, numpy.max([1, numpy.min(rate_series[good4fit])])])
        # ystart = numpy.max([10, numpy.min(rate_series[good4fit])])
        # ax.set_ylim((ystart, 1.3 * numpy.max(rate_series[good4fit])))
        ax.set_ylim((-100, 250))
        ax.plot(read_times, _persistency_plus_signal_fit_fct(bestfit, read_times))

        ax2 = ax.twinx()
        ax2.set_ylabel("linearized read counts")
        ax2.spines['right'].set_color('red')
        # ax2.spines['left'].set_color('blue')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        ax2.scatter(read_times[good4fit], linear_series[good4fit], c='red', s=2)
        ax2.scatter(read_times[~good4fit], linear_series[~good4fit], c='red', facecolors='none', s=2, alpha=0.5)
        ax2.axhline(y=62000, linestyle=":", color='red')

        ax.set_title("S(t) = %.1f + %.1f * exp(-t/%.3f)" % (bestfit[0], bestfit[1], bestfit[2]))

        try:
            plot_fn = "debugplots/debug_y=%04d_x=%04d.png" % (y,x)
            print(plot_fn)
            fig.savefig(plot_fn, bbox_inches='tight')
        except Exception as e:
            print(e)
            pass

        try:
            dmp_fn = "debugplots/dump_y=%04d_x=%04d.png" % (y,x)
            numpy.savetxt(dmp_fn, numpy.array([rate_series, linear_series, uncertainties]).T)
        except:
            pass
        
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

            good_pixel_result = False

            ramp_offset = linearized_cube[0, row, x]

            raw_reads = linearized_cube[:, row, x]
            diff_reads = differential_cube[:, row, x]
            raw_reads = linearized_cube[:, row, x]
            good_data = (raw_reads > 0) & (raw_reads < 62000)
            unsaturated = (raw_reads > 0) & (raw_reads < 62000)
            n_unsaturated = numpy.sum(unsaturated)

            if (full_fit_mask[x]):
                # do a full fit for this pixel
                results = persistency_fit_pixel(
                    differential_cube=differential_cube,
                    linearized_cube=linearized_cube,
                    read_times=read_times,
                    good_data=good_data,
                    x=x, y=row,
                    write_test_plot=write_test_plots,
                )
                best_fit, fit_uncertainties, good4fit = results
                if (best_fit is not None):
                    linebuffer[0:n_pars, x] = best_fit
                    linebuffer[n_pars:2*n_pars, x] = fit_uncertainties

                    integrated_persistency = \
                        best_fit[1] * best_fit[2] * (
                        numpy.exp(-read_times[1]/best_fit[2]) - numpy.exp(-numpy.nanmax(read_times)/best_fit[2]))
                    linebuffer[-3, x] = integrated_persistency
                    good_pixel_result = True

                    # calculate both persistency signal and true signal before saturation
                    p_signal = best_fit.copy()
                    p_signal[1:] = 0
                    only_signal = _persistency_plus_signal_fit_fct2(p=p_signal, read_time=read_times)
                    cumulative_signal = numpy.nancumsum(only_signal)
                    p_persistency = best_fit.copy()
                    p_persistency[0] = 0
                    only_persistency = _persistency_plus_signal_fit_fct2(p=p_persistency, read_time=read_times)
                    cumulative_persistency = numpy.nancumsum(only_persistency)
                    full_data = cumulative_signal + cumulative_persistency
                    # unsaturated = raw_reads <   full_data < 35000
                    total_signal = numpy.nansum(only_signal[unsaturated])
                    total_persistency = numpy.nansum(only_persistency[unsaturated])
                    signal_fraction = total_signal / (total_persistency + total_signal)
                    linebuffer[-1, x] = signal_fraction

                linebuffer[-3, x] = numpy.sum(good4fit)
                linebuffer[-2, x] = 0.

            elif (n_unsaturated > 5):
                # no need for a full fit, just calculate a simple slope
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
                linebuffer[n_pars,x] = _sigma
                linebuffer[-3,x] = numpy.sum(good_data)
                linebuffer[-2,x] = 0.
                good_pixel_result = True

                linebuffer[1, x] = numpy.nanmean(diff_reads[good_data])
                linebuffer[2, x] = numpy.nanmean(diff_reads)

            if (not good_pixel_result or n_unsaturated < 5):
                # Fitting didn't work, and not enough data for a simple
                # slope fit either:
                try:
                    min_rate = numpy.nanmin(diff_reads[unsaturated])
                except:
                    min_rate = 0
                linebuffer[0,x] = min_rate
                linebuffer[n_pars,x] = -99
                linebuffer[-3,x] = -numpy.sum(unsaturated)

                # if ((x == 385 and row == 81) or (x == 387 and row == 95) or (x == 484 and row == 56)):
            if (x >380 and x<390 and row > 90 and row < 100):
                print(x,row, full_fit_mask[x], )
                numpy.savetxt("dump___x=%d_y=%d.txt" % (x,row),
                      numpy.array([diff_reads, raw_reads, good_data.astype(int),
                                   numpy.ones_like(diff_reads)*linebuffer[0,x]]).T)

        # end of loop over all pixels in this row

        persistency_fit[:, row, :] = linebuffer


    return


darktype_GOOD = 0
darktype_COLD = 1
darktype_WARM = 2
darktype_HOT = 3

dump_options = [
    'rawstack',
    'zerosubstack',
    'linearstack',
    'diffstack',
    'cleanstack',
    'darkstack',
    'plain',
    'badmask',
    'weighted',
    'noise',
    'all',
    'fordark',
    'ngoodpixels',
]


def fit_pairwise_slope(times, reads, noise, good_reads=None, plot=False, permplot=True, plottitle=None):
    if (good_reads is None):
        good_reads = numpy.isfinite(reads) & numpy.isfinite(noise) & numpy.isfinite(times) & (times >= 0)

    times = times[good_reads]
    reads = reads[good_reads]
    noise = noise[good_reads]

    # time differences between reads
    dt = times.reshape((-1, 1)) - times.reshape((1, -1))
    df = reads.reshape((-1, 1)) - reads.reshape((1, -1))
    d_noise = noise.reshape((-1, 1)) + noise.reshape((1, -1))

    useful_pairs = (dt > 0) & numpy.isfinite(df)
    rates = (df / dt)[useful_pairs]
    noises = d_noise[useful_pairs]
    #     print(rates.shape)

    try:
        good = numpy.isfinite(rates)

        for it in range(3):
            _stats = numpy.nanpercentile(rates[good], [16, 50, 84])
            _med = _stats[1]
            _sigma = 0.5 * (_stats[2] - _stats[0])
            good = good & (rates > _med - 3 * _sigma) & (rates < _med + 3 * _sigma)
            n_good = numpy.sum(good)
            if (n_good < 10):
                break

        weights = 1. / noises
        weighted = numpy.sum((rates * weights)[good]) / numpy.sum(weights[good])
    except Exception as e:
        weighted,_med,_sigma = numpy.NaN, numpy.NaN, numpy.NaN

    if (plot):
        fig, ax = plt.subplots(tight_layout=True)
        if (permplot):
            perm = numpy.random.permutation(rates.shape[0])
        else:
            perm = numpy.arange(rates.shape[0])
        ax.scatter(numpy.arange(rates.shape[0])[perm][good], rates[good], s=0.5, alpha=.2)
        ax.scatter(numpy.arange(rates.shape[0])[perm][~good], rates[~good], s=0.3, c='red')
        ax.axhline(y=_med)
        ax.axhline(y=_med - _sigma, c='grey')
        ax.axhline(y=_med + _sigma, c='grey')
        ax.axhline(y=bestfit[0], c='orange')
        ax.set_ylim((_med - 10 * _sigma, _med + 10 * _sigma))
        ax.set_title("median/sigma: %.3f +/- %.3f    weighted: %.3f" % (_med, _sigma, weighted))
        if (plottitle is not None):
            fig.suptitle(plottitle)

    return dict(
        weighted=weighted,
        median=_med,
        sigma=_sigma
    )


def worker__reference_pixel_correction(
        shmem_cube_raw, shmem_cube_corrected, cube_shape,
        refpixel_mode,
        jobqueue, workername=None,
):
    """

    :param shmem_cube_raw: shared memory block with raw input data
    :param shmem_cube_corrected: shared memory block to hold corrected output file
    :param cube_shape: shape of data cube
    :param refpixel_mode: what algorithm to use for corrections. supported options are given in TODO: XXX
    :param jobqueue: queue to handle workload balancing. By default reference pixels are corrected read by read
    :param workername: name for worker to use during logging
    :return: No return value, work is done when the jobqueue is empty
    """

    logger = logging.getLogger(workername if workername is not None else "RefPixelWorker")
    logger.debug("Starting worker")

    cube_raw = numpy.ndarray(shape=cube_shape, dtype=numpy.float32,
                             buffer=shmem_cube_raw.buf)
    cube_corrected = numpy.ndarray(shape=cube_shape, dtype=numpy.float32,
                             buffer=shmem_cube_corrected.buf)

    while (True):
        try:
            job = jobqueue.get()
            if (job is None):
                # this is the termination signal
                jobqueue.task_done()
                break
        except queue.Empty as e:
            logger.warning("job queue empty (%s)" % (e))
            break

        slice = job
        logger.debug("Starting reference pixel correction for slice/read %d" % (slice))

        t1 = time.time()
        # get correction from data
        data = cube_raw[slice, :, :]
        corr_image = reference_pixels_to_background_correction(
            data=data,
            edge=1,
            mode=refpixel_mode,
            verbose=False, make_plots=False, debug=False)

        # write results to output shared memory
        cube_corrected[slice, :, :] = data - corr_image
        t2 = time.time()
        logger.debug("Correction for read %d done after %.3f seconds" % (slice, t2-t1))

        jobqueue.task_done()

    logger.debug("Shutting down")
    shmem_cube_corrected.close()
    shmem_cube_raw.close()

def worker__nonlinearity_correction(
        shmem_cube_corrected, shmem_corrections, cube_shape, corrections_shape,
        jobqueue, workername=None,
):
    """

    :param shmem_cube_corrected:
    :param shmem_corrections:
    :param cube_shape:
    :param corrections_shape:
    :param jobqueue:
    :param workername:
    :return:
    """

    logger = logging.getLogger(workername if workername is not None else "NonlinCorrectionWorker")
    logger.debug("Starting worker")

    cube_corrected = numpy.ndarray(shape=cube_shape, dtype=numpy.float32,
                             buffer=shmem_cube_corrected.buf)
    poly_corrections = numpy.ndarray(shape=corrections_shape, dtype=numpy.float32,
                             buffer=shmem_corrections.buf)

    # derive polynomial degree from number of factors
    poly_order = poly_corrections.shape[0] - 1

    while (True):
        try:
            job = jobqueue.get()
            if (job is None):
                # this is the termination signal
                jobqueue.task_done()
                break
        except queue.Empty as e:
            logger.warning("job queue empty (%s)" % (e))
            break

        y = job
        logger.debug("Starting nonlinearity correction for row %d" % (y))

        t1 = time.time()
        # get correction from data

        linecube = numpy.array(cube_corrected[:,y,:])
        linecube[linecube > 50000] = numpy.NaN  # mask out all saturated pixels
        zero_offset = numpy.nanmin(linecube, axis=0)
        # print(zero_offset.shape)
        linecube -= zero_offset

        linefactors = poly_corrections[:,y,:]

        outbuf = numpy.zeros_like(linecube)
        for p in range(poly_order):
            # logger.info("poly-order %d (img^%d)" % (p, poly_order-p))
            outbuf += linefactors[p] * numpy.power(linecube, poly_order-p)

        cube_corrected[:,y,:] = outbuf[:,:]

        t2 = time.time()
        logger.debug("Correction for row %d done after %.3f seconds" % (y, t2-t1))

        jobqueue.task_done()

    logger.debug("Shutting down")
    shmem_corrections.close()
    shmem_cube_corrected.close()





class NIRWALS(object):
    """
    NIRWALS class -- handles all NIRWALS data processing, from reading input frames into data cubes
    to generating the final rate images to writing the final data products to disk.

    """
    mask_SATURATED = 0x0001
    mask_LOW_RATE = 0x0002
    mask_BAD_DARK = 0x0004
    mask_NEGATIVE = 0x0008

    RESULT_EXTENSIONS = ["SCI", "MEDIAN", "NOISE", "NPAIRS", "MAX_T_EXP"]
    N_RESULTS = len(RESULT_EXTENSIONS)

    def __init__(self, fn, max_number_files=-1,
                 saturation=None,
                 saturation_level=50000,
                 saturation_fraction=0.25, saturation_percentile=95,
                 use_reference_pixels='none',
                 algorithm='linreg',
                 mask_saturated_pixels=False,
                 nonlinearity=None,
                 n_cores=0,
                 speedy=False,
                 dumps=None,
                 every=None,
                 correct_gain=False,
                 logger_name=None
                 ):
        """
        Initializer for the NIRWALS reduction class, which takes all configuration parameters needed for operation.

        :param fn: filename of any file from the given read sequence. All other filenames are generated on the fly
        based on this filename and the information in its headers.

        :param max_number_files: Optional: Limit the number of reads to load to minimize memory footprint and/or speed
        up processing

        :param saturation: (Optional) read saturation levels on a pixel-by-pixel level from this file.

        :param saturation_level: Saturation level; all pixels with values above this pixel will be masked out as Infs
        during the data loading stage and ignored during all downstream processing.

        :param saturation_fraction: not implemented

        :param saturation_percentile:  not implemented

        :param use_reference_pixels: Method to use for reference pixel correction.

        :param algorithm: Name of algorithm to use for combining the reference pixel and nonlinearity corrected read
        cubes into the desired rate frame.

        :param correct_gain: Apply amplifier-specific gain correction to the data to yield results in electron/second (
        (rather than the default of ADU counts/second)

        :param mask_saturated_pixels:

        :param nonlinearity: name of filename with nonlinearity correction parameters

        :param n_cores: number of parallel CPU cores during some of the data processign steps.

        :param speedy: select only a subset of input reads to speed up processing; no longer recommended

        :param dumps: list of intermediate data products to write to disk during processing.

        :param every: use only for testing; instead of loading all input reads, only read an evenly spaced subset of
        input reads (i.e. read only every N-th read))

        :param logger_name: name for logger.

        """

        self.fn = fn
        self.filelist = []
        self.logger = logging.getLogger("Nirwals" if logger_name is None else logger_name)

        self.logger.info("Current working directory: %s" % (os.getcwd()))

        self.use_reference_pixels = use_reference_pixels
        self.algorithm = algorithm
        self.image_stack_initialized = False
        self.first_read_subtracted = False
        self.first_read = None
        self.first_header = None
        self.saturation_frame = None
        self.gain = None

        self.nonlin_fn = nonlinearity
        self.nonlinearity_cube = None
        self.nonlinearity_polyorder = -1

        self.alloc_persistency = False

        # store values we may/will need during reduction
        self.max_number_files = -1 if max_number_files is None else max_number_files
        self.saturation_level = saturation_level
        self.saturation_fraction = saturation_fraction
        self.saturation_percentile = saturation_percentile
        self.mask_saturated_pixels = mask_saturated_pixels

        self.shmem_cube_raw = None
        self.shmem_cube_linearized = None
        self.shmem_cube_nonlinearity = None
        self.shmem_cube_results = None

        # Keep track of what intermediate processing steps to save/dump
        self.write_dumps = dumps

        self.header_first_read = None
        self.header_last_read = None

        self.provenance = DataProvenance(
            logger=self.logger,
            track_machine_data=True
        )

        if (saturation is not None):
            if (os.path.isfile(saturation)):
                self.saturation_frame = saturation
                self.logger.info("Reading pixel-by-pixe saturation limits from %s" % (self.saturation_frame))
                sat_hdu = pyfits.open(saturation)
                self.saturation_level = sat_hdu[0].data
                self.logger.info("saturation limits shape: %s" % (str(self.saturation_level.shape)))
                self.provenance.add("saturation-level", saturation)
            else:
                try:
                    _sat = float(saturation)
                    self.saturation_level = _sat
                    self.logger.info("Setting custom saturation: %.2f" % (_sat))
                    self.provenance.add("saturation-level", _sat)
                except:
                    mplog.log_exception()
                    self.logger.warn("Unable to handle custom saturation level: %s (%s)" % (saturation, type(saturation)))
                    pass
        else:
            self.provenance.add("saturation-level", self.saturation_level)
            self.logger.info("Using default saturation level: %.2f" % (self.saturation_level))

        self.n_cores = n_cores if (n_cores is not None and n_cores > 0) else multiprocessing.cpu_count()
        self.logger.info("Using %d CPU cores/threads for parallel processing" % (self.n_cores))

        self.correct_gain = correct_gain

        self.logger.debug("Reading exposure setup")
        self.read_exposure_setup()
        self.speedy = speedy
        self.every = every

        self.logger.debug("Retrieving filelist")
        self.get_full_filelist()

        try:
            self.logger.debug("Allocating shared memory")
            self.allocate_shared_memory()
        except FileExistsError as e:
            self.logger.critical(
                "Unable to allocate shared memory -- most likely causes are it is either left from previous aborted run "
                "(check /dev/shm or the like) or you are running multiple runs in parallel. Clean up files and try again")
            raise(e)

    def nonlinearity_valid(self):
        """
        Check if the specified non-linearity correction file and the data within it is in the right format.

        :return: True if data is valid, False otherwise

        """
        if (self.nonlin_fn is None or not os.path.isfile(self.nonlin_fn)):
            return False
        try:
            self.logger.debug("Nonlinearity file (%s) exists, checking for NONLINPOLY extension and correct dimensions" % (self.nonlin_fn))
            hdu = pyfits.open(self.nonlin_fn)
            nonlin_ext = hdu['NONLINPOLY']
            data = nonlin_ext.data
            if (data is not None and data.ndim == 3 and data.shape[1]==self.ny and data.shape[2]==self.nx):
                polyorder = data.shape[0] - 1
                self.logger.info("Found valid non-linearity correction cube, order %d, in %s" % (polyorder, self.nonlin_fn))
                self.nonlinearity_polyorder = polyorder
                return True
        except:
            pass
        return False

    def allocate_shared_memory(self):
        """
        Pre-allocate shared memory for all datacubes to miniminze RAM footprint as much as possible, while
        providing optimal performance for parallel data processing.
        """

        n_groups = numpy.min([self.n_groups, self.max_number_files]) if self.max_number_files > 0 else self.n_groups
        cube_shape = (n_groups, self.ny, self.nx)
        n_pixels_in_cube = n_groups * self.ny * self.nx
        dummy = numpy.array([], dtype=numpy.float32)
        self.logger.info("Assuming cube dimensions of %s for shared memory (%d bytes/pixel, total %.3f GB per cube)" % (
            str(cube_shape), dummy.itemsize, n_pixels_in_cube*dummy.itemsize/2**30))

        self.logger.info("Allocating shared memory: raw cube")
        self.shmem_cube_raw = multiprocessing.shared_memory.SharedMemory(
            name='raw_datacube', create=True,
            size=(dummy.itemsize * n_pixels_in_cube),
        )
        self.cube_raw = numpy.ndarray(shape=cube_shape, dtype=numpy.float32, buffer=self.shmem_cube_raw.buf)

        # logger.info("Copying datacube to shared memory")
        # _raw = numpy.ndarray(shape=data_shape, dtype=numpy.float32, buffer=shmem_raw.buf)
        # _raw[:, :, :] = raw_cube[:, :, :]

        # allocate buffer for corrected data
        self.logger.info("Allocating shared memory: linearized cube")
        self.shmem_cube_linearized = multiprocessing.shared_memory.SharedMemory(
            name='linearized_datacube', create=True,
            size=(dummy.itemsize * n_pixels_in_cube),
        )
        self.cube_linearized = numpy.ndarray(shape=cube_shape, dtype=numpy.float32, buffer=self.shmem_cube_linearized.buf)
        self.cube_linearized[:, :, :] = 0.

        # Add place-holders for the nonlinearity corrections, to be read later
        self.shmem_cube_nonlinearity = None
        self.cube_nonlinearity = None
        if (self.nonlinearity_valid()):
            nonlin_shape = (self.nonlinearity_polyorder + 1, self.ny, self.nx)
            n_pixels_nonlinearity = (self.nonlinearity_polyorder + 1) * self.ny * self.nx
            self.logger.info("Allocating shared memory: linearized cube")
            self.shmem_cube_nonlinearity = multiprocessing.shared_memory.SharedMemory(
                name='nonlinearity_corrections', create=True,
                size=(dummy.itemsize * n_pixels_nonlinearity),
            )
            self.cube_nonlinearity = numpy.ndarray(shape=nonlin_shape, dtype=numpy.float32, buffer=self.shmem_cube_nonlinearity.buf)
            self.cube_nonlinearity[:, :, :] = 0.
            self.logger.debug("shared meory for nonlinearity corrections initialized")

        self.n_results_dimension = self.N_RESULTS
        self.logger.info("Allocating shared memory for results (ndim=%d)" % (self.n_results_dimension))
        n_pixels_results_cube = self.n_results_dimension * self.nx * self.ny
        self.shmem_cube_results = multiprocessing.shared_memory.SharedMemory(
            name='results_datacube', create=True,
            size=(dummy.itemsize * n_pixels_results_cube),
        )
        self.cube_results = numpy.ndarray(
            shape=(self.n_results_dimension, self.ny, self.nx), dtype=numpy.float32, buffer=self.shmem_cube_results.buf
        )

    def read_exposure_setup(self):
        if (self.fn is None):
            self.logger.critical("Unable to get exposure setup without valid input filename")

        # read the input file as reference file
        self.logger.debug("Reading setup from %s" % (os.path.abspath(self.fn)))
        try:
            self.ref_hdulist = pyfits.open(self.fn)
            self.ref_header = self.ref_hdulist[0].header
        except:
            self.logger.critical("Unable to open input file (%s)" % (os.path.abspath(fn)))
        self.provenance.add("ref-header", os.path.abspath(self.fn))

        # image dimensions
        self.nx = self.ref_header['XSTOP'] - self.ref_header['XSTART'] + 1
        self.ny = self.ref_header['YSTOP'] - self.ref_header['YSTART'] + 1

        # readout settings -- FOR NOW NO MULTIPLE READS/RAMPS SUPPORTED IN A SINGLE SEQUENCE
        self.n_groups = self.ref_header['NGROUPS']
        self.n_ramps = 1 #self.ref_header['NRAMPS']
        self.n_reads = 1 #self.ref_header['NREADS']
        self.n_outputs = self.ref_header['NOUTPUTS']

        # exposure and other times
        self.exptime = self.ref_header['USEREXP'] / 1000.
        self.diff_exptime = self.exptime / self.n_groups

        if (self.correct_gain):
            self.gain = nirwals.data.NirwalsGain(self.ref_header)


    def get_full_filelist(self):
        # get basedir
        fullpath = os.path.abspath(self.fn)
        self.logger.info("absolute path of input file: %s" % (fullpath))

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
        """
        Not implemented; Add additional file during execution for on-the-fly reprocessing.
        :param filename:
        :return:
        """
        return

    def load_all_files(self, max_number_files=None, mask_saturated_pixels=True):
        """
        Load all input files and perform initial data masking.

        :param max_number_files:  Maximum number of files to read to speed up processing and reduce memory footprint.
        Defaults to all files
        :param mask_saturated_pixels: mask pixels above the saturation limit as Infs to exclude them from processing
        down the line; activated by default
        :return: Nothing

        """

        if (max_number_files is None):
            max_number_files = self.max_number_files
        if (self.image_stack_initialized):
            self.logger.debug("stack already initialized, skipping repeat try")
            return

        self._image_stack = []

        # open all frames
        _filelist = self.filelist

        if (self.every is not None):
            # special test mode
            _filelist = _filelist[::self.every]
            self.n_groups = len(_filelist)

        # if (max_number_files > 0):
        #     print("Limiting filelist to %d files" % (max_number_files))
        #     _filelist = _filelist[:max_number_files]

        # setup the data-cube to hold all the data
        if (max_number_files is not None and max_number_files > 0 and self.n_groups > max_number_files):
            self.n_groups = max_number_files
            self.logger.info("Limiting input data to %d read-groups" % (max_number_files))

        # self.image_stack_raw = numpy.full(
        #     (self.n_reads, self.n_groups, self.ny, self.nx),
        #     fill_value=numpy.NaN, dtype=numpy.float32)
        self.raw_read_times = numpy.full((self.n_groups), fill_value=numpy.NaN)

        self.logger.debug("raw image cube dimensions: %s" % (str(self.cube_raw.shape)))


        # TODO: Add proper handling for combined Fowler and up-the-ramp sampling
        for _if, fn in enumerate(_filelist):
            try:
                hdulist = pyfits.open(fn)
                hdr = hdulist[0].header
                imgdata = hdulist[0].data.astype(numpy.float32)
            except Exception as e:
                self.logger.error("Unable to open %s (%s)" % (fn, str(e)))
                continue

            if (self.header_first_read is None):
                self.header_first_read = hdr
            self.header_last_read = hdr

            # hdulist.info()

            img_group = hdr['GROUP'] if self.every is None else _if #hdr['GROUP']//self.every
            img_read = hdr['READ']
            img_exptime = hdr['ACTEXP'] / 1000000. # convert time from raw microseconds to seconds
            self.logger.debug("FN=%s // grp=%d rd=%d exptime=%.4f" % (fn, img_group, img_read, img_exptime))

            if (max_number_files > 0 and img_group > max_number_files):
                self.logger.debug("img-group > max-number-file --> skipping this file")
                continue

            # track input files for provenance
            self.provenance.add("input", fn)

            # mask all saturated pixels if requested
            if (mask_saturated_pixels or True):
                saturation_mask = (imgdata > self.saturation_level) & \
                                  numpy.isfinite(self.saturation_level)
                # no matter what, NEVER mask out reference pixels, otherwise bad things may happen
                saturation_mask[  :  ,   :4] = False # left
                saturation_mask[  :  , -4: ] = False # right
                saturation_mask[  :4 ,   : ] = False # top
                saturation_mask[-4:  ,   : ] = False # bottom
                self.logger.debug("masking out %d saturated pixels" % (
                    numpy.sum(saturation_mask)))
                imgdata[saturation_mask] = numpy.Inf

            self.raw_read_times[img_group-1] = img_exptime
            # self.logger.debug("raw read times: %s" % (str(self.raw_read_times)))

            # self.image_stack_raw[img_read-1, img_group-1, :, :] = imgdata

            # self.logger.debug("Reading read %d into slice %d" % (img_group, img_group-1))
            self.cube_raw[img_group-1, :, :] = imgdata


            # self._image_stack.append(imgdata)
            # if  (self.first_header is None):
            #     self.store_header_info(hdulist[0].header)

            # break

        # calculate the initial image stack
        self.logger.info("#groups=%d // #ramps=%d // #reads=%d" % (self.n_groups, self.n_ramps, self.n_reads))
        if (self.n_groups == 1 and self.n_reads >= 1):
            # this is a plain fowler mode, so calculate pair-wise differences
            self.logger.critical("FOWLER MODE NOT SUPPORTED, NEED N_GROUPS>1 AND N_READS==1")
            # self.logger.info("Using fowler-sampling strategy")
            # self.image_stack = numpy.diff(self.image_stack_raw, axis=0)
            # print("@@@@@@@@@@@@@", self.image_stack.shape)
            # self.read_times = numpy.nanmean(self.raw_read_times, axis=0)
            return

        elif (self.n_groups > 1):
            # up-the-ramp sampling.
            self.logger.info("Using up-the-ramp strategy")
            # self.image_stack = numpy.nanmean(self.image_stack_raw, axis=0)
            # print(self.raw_read_times)
            self.read_times = self.raw_read_times #numpy.nanmean(self.raw_read_times, axis=0)

        else:
            self.logger.critical("No idea what's going here and what to do with this data - HELP!!!!")
            return

        # self.logger.debug("stack before/after: %s --> %s" % (str(self.image_stack_raw.shape), str(self.image_stack.shape)))
        self.logger.debug("read-times: %s" % (str(self.read_times)))
        self.logger.debug("stack shape: %s" % (str(self.cube_raw.shape)))

        #
        # Now also load the nonlinearity corrections from file into shared memory
        #
        if (self.nonlinearity_valid()):
            self.logger.info("Loading nonlinearity correction factors from %s" % (self.nonlin_fn))
            hdu = pyfits.open(self.nonlin_fn)
            nonlin = hdu['NONLINPOLY'].data
            self.cube_nonlinearity[:,:,:] = nonlin[:,:,:]


        # delete raw stack to clean up memory
        # del self.image_stack_raw

        self.image_stack_initialized = True
    def apply_dark_correction(self, dark_fn, dark_mode="differential"):
        """
        Apply dark correction to the input data, as given in counts/second; no longer functional, do not use!

        :param dark_fn:
        :param dark_mode:
        :return:
        """
        self.dark_cube = numpy.zeros_like(self.linearized_cube)
        if (dark_fn is None):
            self.logger.warning("No dark correction requested, skipping")
            return
        elif (not os.path.isfile(dark_fn)):
            self.logger.warning("Dark requested (%s) but not found" % (dark_fn))
            return
        else:
            try:
                # load dark-rate image
                self.logger.info("Loading dark-corrections from %s" % (dark_fn))
                dark_hdu = pyfits.open(dark_fn)
                dark = dark_hdu[1].data
                self.provenance.add("dark", dark_fn)
            except Exception as e:
                self.logger.error("error during dark subtraction:\n%s" % (str(e)))
                return

        if (dark_mode == "linearized"):
            # TODO: CHECK THAT THIS IS CORRECT
            # perform a dark subtraction;
            # dark-current = rate [in cts/sec] * frame-# * exposure-time per frame [in sec]
            self.dark_cube = (numpy.arange(linearized.shape[0], dtype=float).reshape((-1, 1, 1)) + 1) \
                        * self.diff_exptime \
                        * dark.reshape((1, dark.shape[0], dark.shape[1]))
            self.logger.debug("shape of dark cube: %s" % (self.dark_cube.shape))
            self.linearized_cube -= self.dark_cube

        elif (dark_mode == "differential"):
            # if (dark.shape[0] == self.differential_cube.shape[0]):
            #     # dimensions match, nothing to do
            #     dark_correction = dark
            # elif (dark.shape[0] > self.differential_cube.shape[0]):
            #     # we have more dark data then is needed
            #     dark_correction = dark[:self.differential_cube.shape[0]]
            # else:
            #     # we don't have enough data - that's bad, so let's try to extrapolate the rest
            #     n_missing = self.differential_cube.shape[0] - dark.shape[0]
            #     dark_correction = numpy.pad(dark, pad_width=((0,n_missing), (0,0), (0,0)), mode='edge')

            # for bad data do not apply any corrections (better safe than sorry)
            dark[~numpy.isfinite(dark)] = 0.

            self.logger.info("Applying differential dark correction (shape: %s / data: %s)" % (
                str(dark.shape), str(self.differential_cube.shape)))
            # pyfits.PrimaryHDU(data=dark_correction).writeto("deleteme__darkcorr.fits", overwrite=True)
            # print(self.differential_cube.shape, dark.shape)
            # pyfits.PrimaryHDU(data=self.differential_cube).writeto("darkcorrect_before.fits", overwrite=True)

            extra_darks = None
            self.logger.info("Dark correction shape: %s" % (str(dark.shape)))
            for read in range(self.differential_cube.shape[0]):
                if (read >= dark.shape[0]):
                    if (extra_darks is None):
                        extra_darks = numpy.nanmean(dark[-5:, :, :], axis=0)
                    self.differential_cube[read, :, :] -= extra_darks
                    self.logger.debug("Dark correction for read %d from dark extrapolation (only %d reads)" % (read, dark.shape[0]))
                else:
                    self.logger.debug("Dark correction for read %d from master dark" % (read))
                    self.differential_cube[read, :, :] -= dark[read, :, :]
            # self.differential_cube[:,:,] = self.differential_cube[:,:,:] - dark_correction[:,:,:] # * 1e-6)
            # pyfits.PrimaryHDU(data=self.differential_cube).writeto("darkcorrect_after.fits", overwrite=True)

        return
    def dump_data(self, data, fn, datatype="?_default_?", extname=None):
        """
        Helper function to write intermediate data products to file, adding the appropriate header values

        :param data: data frame/cube to write
        :param fn:  output filename
        :param datatype: used to identify type of data for logging
        :param extname: name of extension in output file
        :return:  None

        """
        self.logger.info("Writing %s to %s" % (datatype, fn))

        # all output gets reference header information
        _ext = [pyfits.PrimaryHDU(header=self.ref_header)]

        # add one or all image datasets
        if (type(data) == list):
            for (extname, d) in data:
                _ext.append(pyfits.ImageHDU(data=d, name=extname))
        else:
            _ext.append(pyfits.ImageHDU(data=data, name=extname))

        # add exposure time list

        # add data provenance, just in case

        # now write results
        hdulist = pyfits.HDUList(_ext)
        hdulist.writeto(fn, overwrite=True)

        return
    def apply_reference_pixel_corrections(self):
        """
        compute and derive the reference pixel correction for the raw input data cube. The actual correction is done
        in XXX, this function handles the parallel execution.

        :return: None
        """
        n_reads = self.cube_raw.shape[0]
        self.provenance.add('reference-pixel-mode', self.use_reference_pixels)

        # prepare jobs for each worker
        self.logger.info("Preparing jobs for ref pixel correction workers")
        slice_queue = multiprocessing.JoinableQueue()
        for slice in range(n_reads):
            slice_queue.put(slice)

        # start workers (as many as requested, but not more than jobs available)
        self.logger.info("Starting ref pixel workers")
        n_workers = numpy.min([self.n_cores, n_reads])
        worker_processes = []
        for n in range(n_workers):
            p = multiprocessing.Process(
                target=worker__reference_pixel_correction,
                kwargs=dict(
                    shmem_cube_raw=self.shmem_cube_raw,
                    shmem_cube_corrected=self.shmem_cube_linearized,
                    cube_shape=self.cube_raw.shape,
                    refpixel_mode=self.use_reference_pixels,
                    jobqueue=slice_queue,
                    workername='RefPixelWorker_%03d' % (n+1)
                ),
                daemon=True
            )
            slice_queue.put(None)
            p.start()
            worker_processes.append(p)

        # wait for work to be completed
        self.logger.info("Waiting for ref pixel work to be done")
        slice_queue.join()

        # make sure all workers are shut down
        for p in worker_processes:
            p.join()

        self.logger.info("All ref pixel correction work complete")

        pass

    def reduce(self, dark_fn=None, mask_bad_data=None, mask_saturated_pixels=False, group_cutoff=None,
               algorithm=None):
        """
        Main function called for the instrumental detrending and read fitting.

        :param dark_fn:
        :param mask_bad_data:
        :param mask_saturated_pixels:
        :param group_cutoff:
        :param algorithm: Which algorithm to use to combine multiple reads/groups into the final rate image.

        :return:
        """
        if (algorithm is not None):
            self.algorithm = algorithm

        self.load_all_files(mask_saturated_pixels=mask_saturated_pixels)
        self.logger.info("Done loading all files")

        self.logger.info("Typical interval between reads: %.5f seconds" % (
            numpy.nanmean(numpy.diff(self.read_times[1:-1]))))

        # pyfits.PrimaryHDU(data=self.image_stack).writeto("raw_stack_dmp.fits", overwrite=True)

        # self.reference_corrections_cube = numpy.full_like(self.image_stack, fill_value=0.)

        # reset_frame_subtracted = self.image_stack.copy()
        # if (self.use_reference_pixels != 'none'):
        # pyfits.PrimaryHDU(data=self.cube_raw).writeto("cube_raw.fits", overwrite=True)
        self.dump_save(imgtype='raw')

        self.logger.info("Applying reference pixel corrections [%s]" % (self.use_reference_pixels))
        self.apply_reference_pixel_corrections()
        self.dump_save(imgtype='refpixcorr')

        if (not self.nonlinearity_valid()):
            self.logger.warning("No valid non-linearity correction selected or found, this is not good")
        else:
            self.apply_nonlinearity_corrections()
        self.dump_save(imgtype='linearized')

        if (self.correct_gain):
            self.apply_gain_correction()
            self.provenance.add('gain', self.gain.get_name())
        else:
            self.provenance.add('gain', 'no-gain-correction')


        # self.logger.info("Dumping corrected datacube to file")
        # pyfits.PrimaryHDU(data=self.cube_linearized).writeto("dump_cube.fits", overwrite=True)

        # self.fit_pairwise_slopes(algorithm="rauscher2007")
        self.provenance.add("urg-algorithm", self.algorithm)
        self.fit_pairwise_slopes(algorithm=self.algorithm, group_cutoff=group_cutoff)
        # self.fit_pairwise_slopes(algorithm="pairwise_slopes")

    def dump_save(self, imgtype=None):
        """
        Helper function to write intermediate data products to file. This function handles the correct data
        selection and filename generation, actual writing-to-disk is performed by calling :dump_data:

        :param imgtype: what data product to write

        :return:
        """
        if (not self.write_dumps):
            # no dumps need to be written
            return
        elif ('all' in self.write_dumps):
            pass
        elif (not imgtype in self.write_dumps):
            # this dump is not selected
            return

        bn = self.filebase + "__"
        self.logger.debug("Dump-file basename: %s" % (bn))

        if (imgtype == 'raw'):
            self.dump_data(self.cube_raw, bn+"stack_raw.fits", "RAW")
        elif (imgtype == 'refpixcorr'):
            self.dump_data(self.cube_linearized, bn + "stack_refpixcorr.fits", "REFPIXCORR")
        elif (imgtype == 'linearized'):
            self.dump_data(self.cube_linearized, bn + "stack_linearized.fits", "LINEARIZED")
        else:
            logger.info("Dumping of %s not yet implemented" % (imgtype))

        return

        # for frame_id in range(self.image_stack.shape[0]):
        #     reference_pixel_correction = reference_pixels_to_background_correction(
        #         self.image_stack[frame_id], debug=False, mode=self.use_reference_pixels,
        #     )
        #     self.reference_corrections_cube[frame_id] = reference_pixel_correction
        #     reset_frame_subtracted[frame_id] -= reference_pixel_correction

        # else:
        #     self.logger.info("Subtracting first read from stack")
        #     # self.subtract_first_read()
        #     # apply first-read subtraction
        #     self.reset_frame = self.image_stack[0]
        #     self.reference_corrections_cube[:] = self.reset_frame.reshape((-1, self.ny, self.nx))
        #     # reset_frame_subtracted = self.image_stack - self.reset_frame

        # print(self.read_times)
        self.logger.info("Typical interval between reads: %.3f seconds" % (
            numpy.nanmean(numpy.diff(self.read_times[1:-1]))))

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

        # also initialize a 16-bit mask frame to hold some more info about the data results
        self.logger.debug("Allocating memory for the output mask frame")
        mask_dtype = numpy.int16
        _mask = numpy.zeros((self.ny, self.nx), dtype=mask_dtype)
        self.shmem_image_mask = multiprocessing.shared_memory.SharedMemory(
            create=True, size=_mask.nbytes
        )
        self.image_mask = numpy.ndarray(
            shape=(self.ny, self.nx), dtype=mask_dtype,
            buffer=self.shmem_image_mask.buf
        )
        self.image_mask[:,:] = 0

        # TODO: Add flexibility to have darks either in linear or differential mode
        # for now let's simplify things and assume differential mode only

        # # allocate shared memory for the differential stack and calculate from
        # # the linearized cube
        # # calculate differential stack
        # self.logger.debug("Allocating shared memory for differential cube")
        # self.shmem_differential_cube = multiprocessing.shared_memory.SharedMemory(
        #     create=True, size=self.linearized_cube.nbytes
        # )
        # self.differential_cube = numpy.ndarray(
        #     shape=self.linearized_cube.shape, dtype=numpy.float32,
        #     buffer=self.shmem_differential_cube.buf,
        # )
        # self.logger.debug("differential cube allocated")
        #
        # self.logger.debug("calculating differential cube")
        # self.differential_read_times = numpy.pad(numpy.diff(self.read_times), (1,0))
        # self.differential_cube[:, :, :] = numpy.pad(
        #     numpy.diff(linearized, axis=0), ((1,0),(0,0),(0,0))
        # ) / self.differential_read_times.reshape((-1,1,1))
        # self.logger.debug("diff stack: %s" % (str(self.differential_cube.shape)))
        #
        # self.logger.info("Next up (maybe): dark correction")
        # self.apply_dark_correction(dark_fn)
        #
        # # mask out all saturated and/or otherwise bad samples
        # max_count_rates = -1000 # TODO: FIX THIS numpy.nanpercentile(self.differential_stack, q=self.saturation_percentile, axis=0)
        # # print("max counrates:", max_count_rates.shape)
        #
        # # TODO: implement full iterative outlier rejection here
        # if (mask_bad_data is None):
        #     mask_bad_data = self.mask_BAD_DARK | self.mask_SATURATED | self.mask_LOW_RATE | self.mask_NEGATIVE
        # self.logger.info("Identifying bad/dead/saturated/negative pixels (0x%02x)" % (mask_bad_data))
        # bad_data = numpy.zeros_like(self.image_stack, dtype=bool)
        # if (mask_bad_data is not None and (mask_bad_data & self.mask_SATURATED) > 0):
        #     bad_data = bad_data | (self.image_stack > self.saturation_level)
        # if (mask_bad_data is not None and (mask_bad_data & self.mask_LOW_RATE) > 0):
        #     bad_data = bad_data | (self.differential_cube < self.saturation_fraction * max_count_rates)
        # # if (mask_bad_data is not None and (mask_bad_data & self.mask_BAD_DARK) > 0):
        # #     bad_data = bad_data | (self.dark_cube >= linearized)
        # # if (mask_bad_data is not None and (mask_bad_data & self.mask_NEGATIVE) > 0):
        #     bad_data = bad_data | (linearized < 0)
        #
        # self.bad_data_mask = bad_data
        #
        # # bad_data = (self.image_stack > self.saturation_level) | \
        # #            (self.differential_stack < self.saturation_fraction*max_count_rates) | \
        # #            (dark_cube >= linearized) | \
        # #            (linearized < 0)
        #
        # self.logger.info("Cleaning image cube")
        # self.clean_stack = self.differential_cube.copy()
        # self.clean_stack[bad_data] = numpy.NaN
        # self.clean_stack[0, :, :] = numpy.NaN # mask out the first slice, which is just padding
        # # pyfits.PrimaryHDU(data=self.clean_stack).writeto("darkcorrect___cleanstack.fits", overwrite=True)
        #
        # # calculate a average countrate image
        # self.logger.info("calculating final image from stack")
        # # image7 = numpy.nanmean(self.clean_stack[:7], axis=0)
        # self.reduced_image_plain = numpy.nanmean(self.clean_stack, axis=0)
        # noise = numpy.sqrt(self.image_stack)
        # # pyfits.PrimaryHDU(data=noise).writeto("darkcorrect___noise.fits", overwrite=True)
        # noise[bad_data] = numpy.NaN
        # noise[0, :, :] = numpy.NaN
        # self.inv_noise = numpy.nansum(1./noise, axis=0)
        # self.weighted_mean = numpy.nanmean(self.clean_stack, axis=0) # numpy.nansum(self.clean_stack / noise, axis=0) / self.inv_noise
        # self.noise_image = 1. / self.inv_noise
        # # print(image.shape)

        # self.weighted_mean = numpy.zeros((2048,2048))
        # self.noise_image = numpy.zeros((2048,2048))
        # self.median_image = numpy.zeros((2048,2048))
        # readnoise = 20
        # for y in range(2048):
        #     self.logger.info("Reconstructing image,y=%d" % (y))
        #     for x in range(2048):
        #         raw_reads = self.linearized_cube[:,y,x]
        #         result = fit_pairwise_slope(
        #             times=self.read_times,
        #             reads=raw_reads,
        #             noise=numpy.sqrt(numpy.fabs(raw_reads) + readnoise**2),
        #             good_reads=None,
        #             plot=False
        #         )
        #         self.weighted_mean[y,x] = result['weighted']
        #         self.noise_image[y,x] = result['sigma']
        #         self.median_image[y,x] = result['median']
        #
        # pyfits.PrimaryHDU(data=self.weighted_mean).writeto("safety__weightedmean.fits", overwrite=True)
        # pyfits.PrimaryHDU(data=self.noise_image).writeto("safety__sigma.fits", overwrite=True)
        # pyfits.PrimaryHDU(data=self.median_image).writeto("safety__median.fits", overwrite=True)
        # ratios = linearized / linearized[3:4, :, :]


        return
    def subtract_first_read(self):
        """
        Subtract the first read from all frames; No longer required as this functionality has been included
        as part of the non-linearity processing.

        :return:
        """
        if (self.first_read_subtracted):
            self.logger.debug("First read already subtracted, skipping")
            return

        if (self.first_read is None):
            self.first_read = self.image_stack[0].copy()

        self.image_stack -= self.first_read
        self.first_read_subtracted = True

    def _nonlinearity_fit_fct(self, p, x):
        """
        Helper function for the nonlinearity fitting; generates the model value based on times and fitting parameters.

        :param p:
        :param x:
        :return:
        """
        y = numpy.zeros(x.shape)
        for i in range(p.shape[0]):
            y += p[i] * x**(i+1)
        return y
    def _nonlinearity_fit_err_fct(self, p, x, y):
        """
        Helper function for the nonlinearity function; using the read times and fitting parameters calculates the
        relative deviations between model and data, relative to the estimated data uncertainty (see chi^2 algorithm).
        This function determines what is minimized during least-squares fitting.

        :param p: input fitting parameters
        :param x: input read times/frame number
        :param y: input data
        :return:

        """
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
        """
        No longer implemented; use functionality in nirwals_fit_nonlinearity instead.

        :param ref_frame_id:
        :param max_linear:
        :param make_plot:
        :return:
        """
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
                bad = (raw_series > 63000) | (diffs < 0.3 * max_diff)  | ~numpy.isfinite(raw_series)

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
        """
        Read the nonlinearity coefficients from the specified input frame

        :param nonlin_fn: filename with nonlinearity coefficients.
        :return:
        """
        self.logger.debug("Reading non-linearity: %s (file exists: %s)" % (
            nonlin_fn, os.path.isfile(nonlin_fn)))
        if (os.path.isfile(nonlin_fn)):
            try:
                self.logger.info("Reading nonlinearity corrections from %s" % (nonlin_fn))
                hdulist = pyfits.open(nonlin_fn)
                # hdulist.info()
                nonlinearity_cube = hdulist['NONLINPOLY'].data
                self.logger.debug("CORR shape: %s" % (str(nonlinearity_cube.shape)))
            except Exception as e:
                self.logger.error(str(e))
                return False

        self.nonlin_fn = nonlin_fn
        self.nonlinearity_cube = nonlinearity_cube

    def apply_nonlinearity_corrections(self):
        """
        Apply the non-linearity corrections

        :return:
        """
        # pyfits.PrimaryHDU(data=self.cube_linearized).writeto("cube_before_nonlin.fits", overwrite=True)

        self.logger.info("Starting nonlinearity correction")
        self.provenance.add("non-linearity", os.path.abspath(self.nonlin_fn))

        t1 = time.time()
        jobqueue = multiprocessing.JoinableQueue()
        for y in range(2048):
            jobqueue.put(y)

        # setup and start worker processes
        worker_processes = []
        n_workers = self.n_cores
        for n in range(n_workers):
            p = multiprocessing.Process(
                target=worker__nonlinearity_correction,
                kwargs=dict(shmem_cube_corrected=self.shmem_cube_linearized,
                            shmem_corrections=self.shmem_cube_nonlinearity,
                            cube_shape=self.cube_raw.shape,
                            corrections_shape=self.cube_nonlinearity.shape,
                            jobqueue=jobqueue,
                            workername="NonLinWorker_%03d" % (n+1),
                ),
                daemon=True
            )
            jobqueue.put(None)
            p.start()
            worker_processes.append(p)

        # wait for work to be done
        self.logger.info("Working for all nonlinearity correction jobs to finish")
        jobqueue.join()
        for p in worker_processes:
            p.join()

        t2 = time.time()

        self.logger.info("Non-linearity correction complete after taking %.3f seconds" % (t2-t1))
        # pyfits.PrimaryHDU(data=self.cube_linearized).writeto("cube_after_nonlin.fits", overwrite=True)
        return

    def apply_gain_correction(self):
        self.logger.info("Applying gain correction: %s", self.gain.get_name())
        self.logger.debug("Gain values used: %s", ",".join("%.4f" % g for g in self.gain.get_gains()))
        gain_strip = self.gain.amp_corrections()
        for read in range(self.cube_linearized.shape[0]):
            self.cube_linearized[read,:,:] *= gain_strip
        self.logger.debug("Done with applying GAIN corrections")
        return

    def fit_pairwise_slopes(self, algorithm='pairwise_slopes', group_cutoff=None):
        """
        Perform the cube fitting, combining multiple reads into a single countrate image. The actual fitting functions
        are implemented elsewhere (in nirwals_urg_algorithms), but this function handles all the parallel processing.

        :param algorithm: what algorithm to use.
        :param group_cutoff: limits the number of read-groups to process.
        :return:
        """
        self.logger.info("Start of up-the-ramp fitting (algorithm: %s)" % (algorithm))
        t1 = time.time()
        jobqueue = multiprocessing.JoinableQueue()
        for y in range(4, 2048-4):
            jobqueue.put(y)

        # setup and start worker processes
        worker_processes = []
        n_workers = self.n_cores
        recombine = (algorithm.lower().endswith("+recombine"))
        if (recombine):
            self.logger.info("Activating URG fixing via sequence recombination")
        # else:
        #     self.logger.info("No sequence recombination requested")
        if (algorithm == 'pairwise_slopes'):
            worker_method = worker__fit_pairwise_slopes
        elif (algorithm.startswith('rauscher2007')):
            worker_method = worker__fit_rauscher2007
        elif (algorithm.startswith('linreg')):
            worker_method = worker__fit_linear_regression
        else:
            self.logger.critical("URG algorithm (%s) not implemented")
            return False

        if (group_cutoff is not None):
            self.logger.info("Restricting data analysis to %d groups" % (group_cutoff))

        for n in range(n_workers):
            p = multiprocessing.Process(
                target=worker_method,
                kwargs=dict(
                    shmem_cube_corrected=self.shmem_cube_linearized,
                    shmem_results=self.shmem_cube_results,
                    cube_shape=self.cube_linearized.shape,
                    results_shape=self.cube_results.shape,
                    jobqueue=jobqueue,
                    workername="PairSlopeWorker_%03d" % (n+1),
                    read_times=self.read_times,
                    speedy=self.speedy,
                    group_cutoff=group_cutoff,
                    recombine=recombine,
                ),
                daemon=True
            )
            jobqueue.put(None)
            p.start()
            worker_processes.append(p)

        # wait for work to be done
        self.logger.info("Waiting for all URG-fitting jobs to finish")
        jobqueue.join()
        for p in worker_processes:
            p.join()

        t2 = time.time()

        self.logger.info("Up-the-ramp slope-fitting complete after taking %.3f seconds" % (t2-t1))
        # pyfits.PrimaryHDU(data=self.cube_linearized).writeto("cube_after_nonlin.fits", overwrite=True)
        return

    def fix_final_headers(self):
        """
        Prepare the final image extension headers before writing results to disk.

        :return:
        """
        # Set some basic headers just to make sure
        self.ref_header['TELESCPE'] = "SALT"
        self.ref_header['INSTRUME'] = "NIRWALS"

        # Delete headers we no longer need, since they are not applicable to the full sequence
        self.logger.debug("Deleting unnecessary headers")
        for key in ['READ', 'ACTEXP']:
            try:
                del self.ref_header[key]
            except:
                pass

        # Add some additional headers
        self.logger.debug("Adding more headers to primary header")
        self.ref_header['EXPTIME'] = (-1., "exposure time [seconds]")
        if (self.header_last_read is not None):
            try:
                self.ref_header['EXPTIME'] = (self.header_last_read['ACTEXP'] / 1.e6)
            except Exception as e:
                self.logger.error("Unable to set final exposure time (%s)" % (str(e)))

        # modify some existing keys to match the sequence
        for key in ['SDST-', 'UTC-', 'TIME-']:
            self.ref_header[key+"OBS"] = "unknown"
            self.ref_header[key+"END"] = "unknown"
            try:
                self.ref_header[key+"OBS"] = self.header_first_read[key+'OBS']
                self.ref_header.insert(key=key+"OBS", card=(key+"END", self.header_last_read[key+'OBS']), after=True)
            except KeyError:
                self.logger.error("Unable to correct FITS header %(key)s-OBS and/or %(key)s-END" % dict(key=key))
                continue

        self.ref_header['PUPSTA'] = -99999.99
        self.ref_header['PUPEND'] = -99999.99
        try:
            self.ref_header['PUPSTA'] = self.header_first_read['PUPSTA']
            self.ref_header['PUPEND'] = self.header_last_read['PUPEND']
        except Exception as e:
            self.logger.error("Error while writing PUPIL start/end headers (%s)" % (str(e)))
        return


    def write_results(self, fn=None, flat4salt=False):
        """
        Write the final data product to disk. This includes the rate image itself, along with additional metadata, and
        the data provenance table.

        :param fn: filename for output file.
        :param flat4salt:
        :return:
        """
        # Add/modify FITS headers for output
        self.fix_final_headers()

        if (fn is None):
            fn = os.path.join(self.basedir, self.filebase) + ".reduced.fits"

        # collect all output results
        _list = [pyfits.PrimaryHDU(header=self.ref_header)]
        if (flat4salt):
            self.logger.info("Writing flattened version of output")
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
            self.logger.info("Writing output in MEF format")
            for i,name in enumerate(self.RESULT_EXTENSIONS):
                _list.append(pyfits.ImageHDU(data=self.cube_results[i], name=name))
            # _list.extend([
            #     pyfits.ImageHDU(data=self.cube_results[0], name="SCI"),
            #     pyfits.ImageHDU(data=self.cube_results[2], name='NOISE'),
            #     pyfits.ImageHDU(data=self.cube_results[1], name='MEDIAN'),
            # ])
            # try:
            #     for i,extname in enumerate(persistency_values):
            #         self.logger.debug("Adding MEF extension: %s" % (extname))
            #         _list.append(
            #             pyfits.ImageHDU(
            #                 data=self.persistency_fit_global[i, :, :],
            #                 name=extname)
            #         )
            # except Exception as e:
            #     self.logger.critical(str(e))
            #     pass

        _list.append(pyfits.ImageHDU(header=self.header_first_read, name="READ_FIRST"))
        _list.append(pyfits.ImageHDU(header=self.header_last_read, name="READ_LAST"))

        self.logger.debug("Adding data provenance")
        _list.append(self.provenance.write_as_hdu())

        self.logger.debug("Adding fiber map data")
        _list.append(nirwals.data.get_fibermap())

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
    def fit_signal_with_persistency(self, n_workers=0, previous_frame=None, write_test_plots=False):
        """
        No longer functional, do not use
        :param n_workers:
        :param previous_frame:
        :param write_test_plots:
        :return:
        """
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
            self.provenance.add("persistency-reference", previous_frame)

            threshold = 58000
            need_full_persistency_fit = prev_img > threshold
            # need_full_persistency_fit[400:450, 1650:1700] = True
            # only do persistency fitting near the middle - to speed up testing
            # need_full_persistency_fit[:, :950] = False
            # need_full_persistency_fit[:, 1050:] = False
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
        """
        no longer functional, do not use
        :param x:
        :param y:
        :param debug:
        :param plot:
        :return:
        """
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
        """
        Used only for debugging, do not use

        :param x:
        :param y:
        :param filebase:
        :param cumulative:
        :param differential:
        :param diff_vs_cum:
        :param show_fits:
        :param show_errors:
        :param show_plot:
        :return:
        """
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

        readout_times = numpy.arange(series.shape[0], dtype=float) * self.diff_exptime
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
        """
        Class destructor; this is where we release and clean up shared memory allocation used for parallel processing

        :return:
        """
        self.logger.info("Releasing shared memory allocations")
        # self.logger.debug("Running destructor and cleaning up shared memory")
        # clean up shared memory

        for shmem in [
            self.shmem_cube_raw,
            self.shmem_cube_linearized,
            self.shmem_cube_results,
            self.shmem_cube_nonlinearity,
            ]:
            if (shmem is None):
                continue
            try:
                shmem.close()
                shmem.unlink()
            except Exception as e:
                self.logger.warning("Error releasing shared memory: %s" % (str(e)))
                pass

        self.logger.debug("All done releasing shared memory")
