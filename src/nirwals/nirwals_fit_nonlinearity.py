#!/usr/bin/env python3
import multiprocessing
import queue

import sys
import os
import multiparlog as mplog
import argparse
import logging
import itertools
import time

# for verification
import astropy.io.fits as pyfits
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from nirwals import NIRWALS, DataProvenance


NONLIN_FLAG_OK = 0x00
NONLIN_FLAG_INSUFFICIENT_DATA = 0x01
NONLIN_FLAG_NEGATIVE = 0x02
NONLIN_FLAG_BADSLOPE = 0x04
NONLIN_FLAG_NEGATIVE_REFSLOPE = 0x08
NONLIN_FLAG_FITERROR = 0x10

N_EXTRAS = 8
# Extras are (noted are indices):
#  0: saturation level
#  1: full well depth raw
#  2: full well depth, after nonlinearity correction
#  3: sigma-clipped 1-sigma relative precision median
#  4: sigma-clipped 1-sigma relative precision scatter
#  5: raw 1-sigma rel precision scatter
#  6: max raw read value
#  7: max refpixel corrected read value


def fit_pixel_nonlinearity(
        reads_refpixelcorr, reads_raw, times,
        poly_order=3, ref_level=10000, saturation_level=55000, min_flux=200,
        logger=None, n_iterations=4, normalize_linear_term=True, pixel_id=None,
        return_full=False,
    ):

    if (logger is None):
        logger = logging.getLogger("FitPixelNonlinearity")

    # subtract off any residual offsets (read for t=0 should be 0)
    reads_offset = numpy.nanmin(reads_refpixelcorr)
    reads_refpixelcorr -= reads_offset
    flags = 0

    # define fallback solution with no correction (i.e. corrected = 1.0 x input + 0)
    fallback_solution = numpy.zeros((poly_order+1))
    fallback_solution[-2] = 1.

    # set some defaults in case things go bad
    slope_reflevel = numpy.NaN
    nonlin_bestfit = fallback_solution
    good4fit = None
    t_exp_reflevel = -1
    fit_iterations = -1

    extras = [numpy.NaN]*N_EXTRAS

    # identify suitable pixels, and fit with polynomial of specified degree
    good4fit = (reads_raw < saturation_level) & (reads_refpixelcorr > min_flux)
    n_good4fit = numpy.sum(good4fit)

    try:
        # first, get an idealized target slope for the actual intensity
        f_masked = reads_refpixelcorr.copy()
        if (numpy.nanmax(reads_refpixelcorr) > ref_level):
            # print(times[reads_refpixelcorr > ref_level])
            t_exp_reflevel = numpy.nanmin(times[reads_refpixelcorr > ref_level])
            f_masked[(f_masked < ref_level) | ~numpy.isfinite(f_masked)] = 1e9
            n_read_reflevel = numpy.argmin(f_masked)  # t_reads[darksub_diffuse > 10000])
        else:
            t_exp_reflevel = numpy.nanmax(times)
            n_read_reflevel = numpy.max(numpy.arange(reads_refpixelcorr.shape[0])[numpy.isfinite(reads_refpixelcorr)])

        slope_reflevel = reads_refpixelcorr[n_read_reflevel] / times[n_read_reflevel]
        logger.debug("time to %d counts @ read %d w/o dark (slope: %.1f)" % (
            t_exp_reflevel, n_read_reflevel, t_exp_reflevel))

        # identify suitable pixels, and fit with polynomial of specified degree
        if (n_good4fit < 10):
            flags |= NONLIN_FLAG_INSUFFICIENT_DATA
            raise ValueError("Insufficient data to fit (only %d samples)" % (n_good4fit))
        elif (n_good4fit > 50):
            # only if we have enough data skip the first few reads -- these
            # might be affected by a reset anomaly and thus are less trustworthy
            good4fit[times < 7.] = False


        nonlin_bestfit = [1.00, 0.00]
        fit_iterations=[]
        for iteration in range(n_iterations):

            if (nonlin_bestfit[-2] < 0):
                flags |= NONLIN_FLAG_NEGATIVE_REFSLOPE
                raise ValueError("Something is seriously amiss, reference slope is negative")

            # apply corrections from prior iteration
            reads_refpixelcorr += nonlin_bestfit[-1]
            slope_reflevel /= nonlin_bestfit[-2]

            nonlin_results = numpy.polyfit(
                x=reads_refpixelcorr[good4fit],
                y=(times * slope_reflevel)[good4fit],
                deg=poly_order,
                full=True
            )
            nonlin_bestfit = nonlin_results[0]
            fit_iterations.append(nonlin_bestfit)

        # check if ANY of the correction values turn negative
        inp = numpy.arange(min_flux, numpy.nanmax(reads_refpixelcorr)) #50000)
        out = numpy.polyval(nonlin_bestfit, inp)
        if ((out < nonlin_bestfit[-1]).any()):
            # numpy.savetxt("nonlin_negatives_%s.txt" % pixel_id.replace(" ",""), numpy.array([inp, out]).T)
            flags |= NONLIN_FLAG_NEGATIVE
            raise ValueError("Illegal solution, corrected value < 0")

        # also check if first derivative is always positive (i.e. corrected counts monotoneously increase with increasing observed counts)
        bestfit_der1 = numpy.polyder(nonlin_bestfit, m=1)
        slopes = numpy.polyval(bestfit_der1, inp)
        if ((slopes <= 0).any()):
            # numpy.savetxt("nonlin_badslope_%s.txt" % pixel_id.replace(" ",""), numpy.array([inp, out, slopes]).T)
            # numpy.savetxt("nonlin_badslope_%s.raw" % pixel_id.replace(" ",""), numpy.array([
            #     times, reads_raw, reads_refpixelcorr, numpy.polyval(nonlin_bestfit, reads_refpixelcorr)]).T)
            flags |= NONLIN_FLAG_BADSLOPE
            raise ValueError("Illegal solution, negative slope detected")

        # print("slope at signal 0:", numpy.polyval(bestfit_der1, 0))
        # print("mean slope 0-1K", numpy.mean(numpy.polyval(bestfit_der1, numpy.linspace(0, 1000, 100))))

        # now convert all poly coefficients such that the linear term is x1.00
        if (normalize_linear_term):
            p1 = nonlin_bestfit[-2]
            for p in range(poly_order):
                nonlin_bestfit[-(p + 2)] /= numpy.power(p1, p + 1)
            # print('corrected:', nonlin_bestfit)

    except ValueError as e:
        logger.debug("Unable to fit non-linearity for %s (%s)" % (pixel_id, str(e)))

    except Exception as e:  # numpy.linalg.LinAlgError as e:
        flags |= NONLIN_FLAG_FITERROR
        logger.debug("Unable to fit non-linearity for %s (%s)" % (pixel_id, str(e)))


    #
    # calculate all extras
    #
    extras[0] = saturation_level
    extras[6] = numpy.nanmax(reads_raw)
    extras[7] = numpy.nanmax(reads_refpixelcorr)

    # 1 sigma precision
    if (flags == 0): # only if all fitting worked out
        try:
            corrected_reads = numpy.polyval(nonlin_bestfit, reads_refpixelcorr)
            simple_slope = (times * slope_reflevel)
            # logger.debug("1sigma precision: %s", corrected_reads.shape, simple_slope.shape)
            relative_error = corrected_reads / simple_slope
            good_rel_error = numpy.isfinite(relative_error) & good4fit

            extras[5] = numpy.nanstd(relative_error[good_rel_error])
            # numpy.savetxt("relerror.dump", relative_error)
            # numpy.savetxt("relerror2.dump", relative_error[good_rel_error])

            for iteration in range(3):
                stats = numpy.nanpercentile(relative_error[good_rel_error], [16,50,84])
                # print("rel slope stats:", stats)
                _med = stats[1]
                _sigma = 0.5*(stats[2]-stats[0])
                good_rel_error = good_rel_error & (relative_error > (_med-3*_sigma)) & (relative_error < (_med+3*_sigma))
            extras[3] = _med
            extras[4] = _sigma
        except Exception as e:
            mplog.report_exception(e)
            pass

    if (flags != NONLIN_FLAG_OK):
        nonlin_bestfit = fallback_solution

    if (not return_full):
        return nonlin_bestfit, slope_reflevel, flags, extras

    return dict(
        bestfit=nonlin_bestfit,
        refslope=slope_reflevel,
        good4fit=good4fit,
        t_exp_reflevel=t_exp_reflevel,
        fit_iterations=fit_iterations,
        flags=flags,
        extras=extras,
    )


def fit_pixel_nonlinearity_autosaturation(
        reads_refpixelcorr, reads_raw, times,
        poly_order=3, ref_level=10000, saturation_level=55000, min_flux=200,
        logger=None, n_iterations=4, normalize_linear_term=True, pixel_id=None,
        return_full=False, optimize_saturation=True, optimize_dt=10, plot_fn=None,
        debug=False,
):

    if (not optimize_saturation):
        # no need to optimize anything, defer this solution back to the simple function
        return fit_pixel_nonlinearity(
            reads_refpixelcorr, reads_raw, times,
            poly_order=poly_order,
            ref_level=ref_level, saturation_level=saturation_level,
            min_flux=min_flux, logger=logger, n_iterations=n_iterations,
            normalize_linear_term=normalize_linear_term, pixel_id=pixel_id,
            return_full=return_full
        )

    if (logger is None):
        logger = logging.getLogger("FitPixelNonlinearity")

    # subtract off any residual offsets (read for t=0 should be 0)
    reads_offset = numpy.nanmin(reads_refpixelcorr)
    reads_refpixelcorr -= reads_offset
    flags = 0

    for iteration in range(3):

        nonlin_bestfit, slope_reflevel, flags, extras = fit_pixel_nonlinearity(
            reads_refpixelcorr=reads_refpixelcorr,
            reads_raw=reads_raw,
            times=times,
            poly_order=poly_order,
            ref_level=ref_level,
            saturation_level=saturation_level,
            pixel_id=pixel_id,
            min_flux=min_flux, logger=logger, n_iterations=n_iterations,
            normalize_linear_term=normalize_linear_term,
            return_full=False,
        )

        corr = numpy.polyval(nonlin_bestfit, reads_refpixelcorr)
        err = numpy.sqrt(reads_refpixelcorr * 2 + 150 ** 2)

        dt = optimize_dt
        diff_times = times[:-dt]
        diff_times_index = numpy.arange(diff_times.shape[0])
        diff_flux = reads_refpixelcorr[dt:] - reads_refpixelcorr[:-dt]
        diff_corr_flux = corr[dt:] - corr[:-dt]

        ref_countflux = slope_reflevel * dt * 0.72775

        # find time when the differential corrected flux drops below the reference flux
        saturated = diff_corr_flux < 0.5 * ref_countflux
        n_saturated_reads = numpy.sum(saturated)
        if (n_saturated_reads <= 0):
            # no saturated reads found
            saturation_time = numpy.max(times)
            saturation_time_index = times.shape[0]-1
        else:
            saturation_time = numpy.min(diff_times[saturated])
            # saturation_time_index = numpy.min(
            #     numpy.array(numpy.arange(diff_times.shape[0]))[diff_times > saturation_time])
            post_saturation = diff_times > saturation_time
            if (numpy.sum(post_saturation) <= 0):
                saturation_time_index = times.shape[0]-1
            else:
                saturation_time_index = numpy.min(diff_times_index[post_saturation])

        logger.debug("saturation time/index: %.1f / %d" % (saturation_time, saturation_time_index))
        new_saturation_level = reads_raw[saturation_time_index]
        if (debug):
            print(iteration, new_saturation_level, n_saturated_reads, saturation_time_index)
            numpy.savetxt("debugdump_diff_%s_it%d.txt" % (debug, iteration), numpy.array([
                diff_times, diff_flux
            ]).T)
            numpy.savetxt("debugdump_%s_it%d.txt" % (debug, iteration), numpy.array([
                times, reads_refpixelcorr
            ]).T)

        full_well_depth_raw = reads_refpixelcorr[saturation_time_index]
        full_well_depth_corrected = corr[saturation_time_index]
        logger.debug("Updated extras: saturation=%.1f    fwd-raw=%.1f   fwd-corr=%.1f",
                     new_saturation_level, full_well_depth_raw, full_well_depth_corrected)

        if (debug):
            plot_fn = "debug_%s_it%d.png" % (debug, iteration)

        if (plot_fn is not None):

            fig, axs = plt.subplots(ncols=2, figsize=(12, 5))
            ax = axs[0]
            ax.scatter(diff_times, diff_corr_flux, s=0.3, label='raw')
            ax.scatter(diff_times[saturated], diff_corr_flux[saturated], s=0.3, c='red')
            print(slope_reflevel)

            ax.axhline(ref_countflux)

            ax = axs[1]
            ax.scatter(times, reads_raw, s=0.3, label='raw', c='blue')
            ax.axhline(saturation_level, ls="--")
            ax.axhline(new_saturation_level)
            ax.axhline(full_well_depth_raw, ls=':', c='green')
            ax.axhline(full_well_depth_corrected, ls='-.', c='orange')
            ax.scatter(times, reads_refpixelcorr, s=0.3, label='reads', c='green')
            print(reads_refpixelcorr[0])
            # print(err)
            ax.plot(times, corr, label='corr', c='orange')
            # ax.plot(times, refpixcorr[:,y,x]-err)
            print(slope_reflevel)

            # ax.set_xlim(-2,30)
            # ax.set_ylim(-100,3000)

            # ax.set_xlim(-2,150)
            # ax.set_ylim(-100,10000)

            ax.set_ylim(0e3, 70e3)

            ax.plot(times, times * slope_reflevel, lw=5, alpha=0.3, label='ref-slope')
            ax.legend()

            fig.suptitle("Iteration %d: %.2f --> %.2f" % (iteration + 1, saturation_level, new_saturation_level))
            if (plot_fn == True):
                fig.show()
            else:
                fig.savefig(plot_fn)
            # plt.close(fig)

        if (n_saturated_reads > 0):
            # nothing saturated found, no need to refine saturation any further
            pass
            #break

        saturation_level = new_saturation_level

    extras[0] = saturation_level
    extras[1] = full_well_depth_raw
    extras[2] = full_well_depth_corrected
    if (debug): print(saturation_level, extras)

    if (not return_full or True):
        return nonlin_bestfit, slope_reflevel, flags, extras

    # TODO: fix this here
    return dict(
        bestfit=nonlin_bestfit,
        refslope=slope_reflevel,
        good4fit=good4fit,
        t_exp_reflevel=t_exp_reflevel,
        fit_iterations=fit_iterations,
        flags=flags,
        extras=extras,
    )


def nonlinfit_worker(jobqueue, resultqueue, times,
                     shmem_cube_raw, shmem_cube_refpixelcorr, shmem_cube_nonlinpoly, cube_shape,
                     shmem_flags, shmem_extras,
                     poly_order=3, ref_level=10000, saturation_level=55000,
                     optimize_saturationlevel=True,
                     workername="NonLinFitWorker"):

    logger = logging.getLogger(workername)
    logger.debug("Starting worker %s" % (workername))

    # make the shared memory available as numpy arrays
    cube_raw = numpy.ndarray(
        shape=cube_shape, dtype=numpy.float32,
        buffer=shmem_cube_raw.buf
    )
    cube_refpixelcorr = numpy.ndarray(
        shape=cube_shape, dtype=numpy.float32,
        buffer=shmem_cube_refpixelcorr.buf
    )
    cube_nonlinpoly = numpy.ndarray(
        shape=(poly_order+1, cube_shape[1], cube_shape[2]), dtype=numpy.float32,
        buffer=shmem_cube_nonlinpoly.buf
    )
    nonlin_flags = numpy.ndarray(
        shape=(cube_shape[1], cube_shape[2]), dtype=numpy.int16,
        buffer=shmem_flags.buf
    )
    nonlin_extras = numpy.ndarray(
        shape=(N_EXTRAS, cube_shape[1], cube_shape[2]), dtype=numpy.float32,
        buffer=shmem_extras.buf
    )


    while(True):
        t1 = time.time()
        try:
            job = jobqueue.get(timeout=1)
        except (queue.Empty, ValueError) as e:
            logger.warning("Timeout error while waiting for jobs")
            job = None

        if (job is None):
            jobqueue.task_done()
            break

        x,y = job
        reads_refpixelcorr = cube_refpixelcorr[:,y,x]
        reads_raw = cube_raw[:,y,x]

        nonlin_bestfit, slope_reflevel, flags, extras = (
            fit_pixel_nonlinearity_autosaturation(
                reads_refpixelcorr=reads_refpixelcorr,
                reads_raw=reads_raw,
                times=times,
                poly_order=poly_order,
                ref_level=ref_level,
                saturation_level=saturation_level,
                logger=logger,
                pixel_id="x,y= %4d,%4d" % (x+1, y+1),
                optimize_saturation=optimize_saturationlevel,
                debug=False, #"x=%d_y=%d" % (x,y) if x==y else False,
            )
        )

        t2 = time.time()

        if (x == y): print(" ".join(["%.2f" % f for f in extras]))

        # store the results in shared memory for later output
        cube_nonlinpoly[:, y, x] = nonlin_bestfit
        nonlin_flags[y,x] = flags
        nonlin_extras[:,y,x] = extras

        #resultqueue.put((x, y, nonlin_bestfit, t2 - t1))
        resultqueue.put((x, y, t2 - t1, slope_reflevel))

        jobqueue.task_done() 

    shmem_cube_raw.close()
    shmem_cube_refpixelcorr.close()
    shmem_cube_nonlinpoly.close()
    logger.debug("Shutting down worker %s" % (workername))



def main():

    mplog.setup_logging(debug_filename="nirwals_debug.log",
                        log_filename="nirwals_fit_nonlinearity.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NirwalsFitNonlinearity")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--saturation", dest="saturation", default=55000,
                         help="saturation value/file")
    cmdline.add_argument("--no_optimize", dest="optimize_saturation", default=True,
                         action='store_false', help="optimize saturation level as part of fitting")
    cmdline.add_argument("--reflevel", dest="reflevel", default=10000, type=float,
                         help="saturation value/file")
    cmdline.add_argument("--ncores", dest="n_cores", type=int, default=multiprocessing.cpu_count(),
                         help="number of CPU cores to use for parallel fitting")
    cmdline.add_argument("--refpixel", dest="ref_pixel_mode", default='blockyslope2',
                         help="reference pixels mode [default: NO]")
    cmdline.add_argument("--verify", dest="verify", default=False, action='store_true',
                         help="verify results rather than fitting coefficients")
    cmdline.add_argument("--pixels", dest="pixels", type=str,
                         help="list of pixel coordinates or file with coordinates")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    fn = args.files[0]
    saturation_fn = args.saturation
    logger.info("Saturation level optimization: %s" % ("ON" if args.optimize_saturation else "OFF"))

    logger.info("Initializing data")
    rss = NIRWALS(fn, saturation=68000, #saturation_fn,
                  max_number_files=args.max_number_files,
                  use_reference_pixels=args.ref_pixel_mode,)
    logger.info("Reading files")
    rss.load_all_files()
    rss.fix_final_headers()

    logger.info("Applying reference pixel corrections")
    rss.apply_reference_pixel_corrections()

    logger.info("Allocating shared memory for output")
    poly_order = 5

    logger.info("Dumping intermediate products for inspection")
    rss.write_dumps = 'all'
    rss.filebase = "fit_nonlinearity"
    # rss.dump_save(imgtype='raw')
    # rss.dump_save(imgtype='refpixcorr')

    dummy = numpy.array([], dtype=numpy.float32)
    dummy_int = numpy.array([], dtype=numpy.int16)
    cube_shape = rss.cube_raw.shape
    shmem_nonlinpoly = multiprocessing.shared_memory.SharedMemory(
        name='nonlinpoly', create=True,
        size=(dummy.itemsize * (poly_order+1) * cube_shape[1] * cube_shape[2]),
    )
    result_nonlinpoly = numpy.ndarray(
        shape=(poly_order+1, cube_shape[1], cube_shape[2]),
        dtype=numpy.float32, buffer=shmem_nonlinpoly.buf)
    # by default apply no correction, just copy the data
    result_nonlinpoly[:,:,:] = 0.
    result_nonlinpoly[-2,:,:] = 1.

    # allocate memory for the flags frame
    shmem_nonlinpoly_flags = multiprocessing.shared_memory.SharedMemory(
        name='nonlinpoly_flags', create=True,
        size=(dummy_int.itemsize * cube_shape[1] * cube_shape[2]),
    )
    result_nonlinpoly_flags = numpy.ndarray(
        shape=(cube_shape[1], cube_shape[2]),
        dtype=numpy.int16, buffer=shmem_nonlinpoly_flags.buf)

    # allocate memory for extra data (stats, precision, saturation limit, etc)
    shmem_nonlinpoly_extras = multiprocessing.shared_memory.SharedMemory(
        name='nonlinpoly_extras', create=True,
        size=(N_EXTRAS * dummy.itemsize * cube_shape[1] * cube_shape[2]),
    )
    result_nonlinpoly_extras = numpy.ndarray(
        shape=(N_EXTRAS, cube_shape[1], cube_shape[2]),
        dtype=numpy.float32, buffer=shmem_nonlinpoly_extras.buf)

    logger.info("Also generating %d extras alongside the nonlinearity coefficients" % (N_EXTRAS))

    if (not args.verify):
        # rss.fit_nonlinearity(ref_frame_id=4, make_plot=False)

        logger.info("Distributing work for parallel processing")
        jobqueue = multiprocessing.JoinableQueue()
        resultqueue = multiprocessing.Queue()
        times = rss.raw_read_times
        shape = rss.cube_raw.shape

        out_processing_time = numpy.full((2048,2048), fill_value=numpy.NaN)
        out_refslope = numpy.full((2048, 2048), fill_value=numpy.NaN)

        # logger.info("Starting to fill queue")
        n_jobs = 0
        refpixels = 4
        for x,y in itertools.product(range(refpixels,2048-refpixels), range(refpixels, 2048-refpixels)):

            # while we use cube_linearized, we did not actually apply any nonlinearity corrections
            jobqueue.put((x, y))
            n_jobs += 1

        logger.info("Done with filling queue")
        logger.debug("STACK: %d" % (rss.cube_raw.shape[0]))

        worker_processes = []
        for n in range(args.n_cores):
            p = multiprocessing.Process(
                target= nonlinfit_worker,
                kwargs=dict(jobqueue=jobqueue,
                            resultqueue=resultqueue,
                            times=rss.raw_read_times[:rss.cube_raw.shape[0]],
                            shmem_cube_raw=rss.shmem_cube_raw,
                            shmem_cube_refpixelcorr=rss.shmem_cube_linearized,
                            shmem_cube_nonlinpoly=shmem_nonlinpoly,
                            shmem_flags=shmem_nonlinpoly_flags,
                            shmem_extras=shmem_nonlinpoly_extras,
                            cube_shape=rss.cube_raw.shape,
                            poly_order=poly_order,
                            ref_level=args.reflevel,
                            saturation_level=args.saturation,
                            optimize_saturationlevel=args.optimize_saturation,
                            workername="Worker_%03d" % (n+1)),
                daemon=True
            )
            jobqueue.put(None)
            p.start()
            worker_processes.append(p)

        # wait for all work to be done
        logger.info("Waiting for parallel fitting to complete")
        jobqueue.join()

        # gather results
        logger.info("Gathering results")
        output_cube = numpy.full((poly_order+1,2048,2048), fill_value=numpy.NaN)
        for n in range(n_jobs):
            (x,y,cpu_time,reflevel) = resultqueue.get()
            out_processing_time[y,x] = cpu_time
            out_refslope[y,x] = reflevel
            # output_cube[:,y,x] = polyfit
            # print(polyfit)

        # make sure all processes are shut down
        for p in worker_processes:
            p.join()

        if (args.nonlinearity_fn is None):
            out_fn = "nonlinpoly.fits"
        else:
            out_fn = args.nonlinearity_fn

        img_hdu = pyfits.ImageHDU(data=result_nonlinpoly, name="NONLINPOLY")
        img_hdu.header['POLYORDR'] = poly_order
        img_hdu.header['REFLEVEL'] = reflevel
        img_hdu.header['REFPIXEL'] = args.ref_pixel_mode,

        flag_hdu = pyfits.ImageHDU(data=result_nonlinpoly_flags, name="FLAGS")
        flag_hdu.header['FLAG_x00'] = "ok"
        flag_hdu.header['FLAG_x01'] = "insufficient data"
        flag_hdu.header['FLAG_x02'] = "correction turns negative"
        flag_hdu.header['FLAG_x04'] = "correction slope < 0"
        flag_hdu.header['FLAG_x08'] = "negative reference slope"
        flag_hdu.header['FLAG_x10'] = "fitting error"
        output_hdulist = pyfits.HDUList([
            pyfits.PrimaryHDU(header=rss.ref_header),
            img_hdu,
            flag_hdu,
            pyfits.ImageHDU(data=result_nonlinpoly_extras[0], name="SATURATION_LEVEL"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[1], name="FULL_WELL_DEPTH_RAW"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[2], name="FULL_WELL_DEPTH_CORRECTED"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[3], name="PRECISION_MEDIAN_CLIPPED"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[4], name="PRECISION_SIGMA_CLIPPED"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[5], name="PRECISION_SIGMA_FULL"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[6], name="MAX_READ_RAW"),
            pyfits.ImageHDU(data=result_nonlinpoly_extras[7], name="MAX_READ_REFPIXCORR"),
            rss.provenance.write_as_hdu()
        ])
        logger.info("Writing correction coefficients to output FITS (%s)" % (out_fn))
        output_hdulist.info()

        output_hdulist.writeto(out_fn, overwrite=True)

        pyfits.PrimaryHDU(data=out_processing_time).writeto(out_fn[:-5]+"__cputime.fits", overwrite=True)
        pyfits.PrimaryHDU(data=out_refslope).writeto(out_fn[:-5]+"__refslope.fits", overwrite=True)
        logger.info("All done!")

    else:
        # rss.reduce(dark_fn=None,)
        # Verify rather than fit results
        coeff_hdu = pyfits.open(args.nonlinearity_fn)
        coeffs = coeff_hdu[0].data
        print(coeffs.shape)

        # Read all pixel coordinates, either from command line or @file
        pixels = []
        if (args.pixels.startswith("@")):
            with open(args.pixels[1:]) as f:
                lines = f.readlines()
                for l in lines:
                    if (l.strip().startswith("#")):
                        continue
                    items = l.split()
                    xy = [int(round(float(x))) for x in items[0:2]]
                    pixels.append(xy)
        else:
            pairs = args.pixels.split(":")
            for p in pairs:
                items = p.split(",")
                xy = [int(round(float(x))) for x in items[0:2]]
                pixels.append(xy)
        print(pixels)

        with PdfPages("nonlinearity_verification.pdf") as pdf:
            for xy in pixels:
                [x,y] = xy
                print(xy,x,y)

                fig, axs = plt.subplots(ncols=2, figsize=(8,4), tight_layout=True)
                # axs = fig.add_subplot(col)

                raw_sequence = rss.cube_raw[:,y,x]
                raw0 = raw_sequence - numpy.nanmin(raw_sequence)
                read_number = numpy.arange(raw_sequence.shape[0])
                times = rss.raw_read_times[:]


                # find closest point to 5000 counts
                closest = numpy.argmin(numpy.fabs(raw0-10000))
                ref_countrate = raw0[closest] / times[closest]
                ref_counts = times * ref_countrate

                poly = coeffs[:, y,x]

                corrected = numpy.polyval(poly, raw0)
                # print(raw_sequence.shape)
                axs[0].plot(times, ref_counts, alpha=0.3, linewidth=3, color="#808080")
                axs[0].scatter(times, raw0, s=2, alpha=0.2, c='blue', label='raw')
                axs[0].plot(times, raw0, 'b-', linewidth=1, c='blue')
                axs[0].scatter(times, corrected, c='orange', label='linearized', s=1)
                axs[0].legend(loc='upper left')

                maxy = numpy.min([numpy.nanmax(raw0), 65000])
                maxt = numpy.nanmax(times)
                axs[0].set_ylim((-0.03*maxy,1.04*maxy))
                axs[0].set_xlim((-0.03*maxt,1.03*maxt))
                axs[0].set_xlabel("Integration time [seconds]")
                axs[0].set_ylabel("net signal [counts]")


                nl = (raw0-ref_counts)/ref_counts
                max_nl = numpy.nanmax([-0.3, numpy.nanmin(nl)])
                print(max_nl)
                axs[1].scatter(raw0/1e3, nl, s=2, alpha=0.2, c='blue', label='raw')
                axs[1].plot(raw0/1e3, nl, 'b-', linewidth=1, c='blue')
                axs[1].scatter(raw0/1e3, (corrected-ref_counts)/ref_counts, c='orange', label='linearized', s=1)

                axs[1].grid(alpha=0.2)
                axs[1].legend(loc='upper right')
                axs[1].set_ylim((max_nl,0.05)) #(-0.03*maxy,1.04*maxy))
                axs[1].set_xlim((0, 1.03*maxy/1e3)) #((-0.03*maxt,1.03*maxt))
                axs[1].set_xlabel("Raw counts [x1000 counts]")
                axs[1].set_ylabel("non-linearity [(raw-corrected)/corrected]")



                fig.suptitle("Pixel %d , %d" % (x,y))

                pdf.savefig(fig)
                # fig.savefig()
                #

    # Release shared memory
    del rss
    for shmem in [shmem_nonlinpoly, shmem_nonlinpoly_flags, shmem_nonlinpoly_extras]:
        shmem.close()
        shmem.unlink()



if __name__ == "__main__":
    main()