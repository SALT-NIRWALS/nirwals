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

from nirwals import NIRWALS



def nonlinfit_worker(jobqueue, resultqueue, times,
                     shmem_cube_raw, shmem_cube_refpixelcorr, shmem_cube_nonlinpoly, cube_shape,
                     poly_order=3, ref_level=10000, saturation_level=55000, workername="NonLinFitWorker"):

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

        # subtract off any residual offsets (read for t=0 should be 0)
        reads_offset = numpy.nanmin(reads_refpixelcorr)
        reads_refpixelcorr -= reads_offset

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
            logger.debug("time to %d counts @ read %d w/o dark (slope: %.1f)", t_exp_reflevel, n_read_reflevel,
                  t_exp_reflevel)

            # identify suitable pixels, and fit with polynomial of specified degree
            good4fit = reads_raw < saturation_level
            if (numpy.sum(good4fit) > 50):
                # only if we have enough data skip the first few reads -- these
                # might be a affected by a reset anomaly and thus are less trustworthy
                good4fit[times < 7.] = False


            nonlin_bestfit = [1.00, 0.00]
            for iteration in range(4):
                # apply corrections from prior iteration
                reads_refpixelcorr += nonlin_bestfit[-1]
                slope_reflevel /= nonlin_bestfit[-2]

                nonlin_results = numpy.polyfit(reads_refpixelcorr[good4fit], (times*slope_reflevel)[good4fit],
                                               deg=poly_order, full=True)
                nonlin_bestfit = nonlin_results[0]

            # now convert all poly coefficients such that the linear term is x1.00
            p1 = nonlin_bestfit[-2]
            for p in range(poly_order):
                nonlin_bestfit[-(p+2)] /= numpy.power(p1, p+1)
            # print('corrected:', nonlin_bestfit)

        except Exception as e: # numpy.linalg.LinAlgError as e:
            logger.warning("Unable to fit non-linearity for x=%d  y=%d (%s)" % (x,y, str(e)))
            nonlin_bestfit = numpy.full((poly_order+1), fill_value=numpy.NaN)


        t2 = time.time()

        cube_nonlinpoly[:, y, x] = nonlin_bestfit
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
    cmdline.add_argument("--reflevel", dest="reflevel", default=10000, type=float,
                         help="saturation value/file")
    cmdline.add_argument("--ncores", dest="n_cores", default=multiprocessing.cpu_count(),
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
    logger.info("Initializing data")
    rss = NIRWALS(fn, saturation=saturation_fn,
                  max_number_files=args.max_number_files,
                  use_reference_pixels=args.ref_pixel_mode,)
    logger.info("Reading files")
    rss.load_all_files()

    logger.info("Applying reference pixel corrections")
    rss.apply_reference_pixel_corrections()

    logger.info("Allocating shared memory for output")
    poly_order = 5

    dummy = numpy.array([], dtype=numpy.float32)
    cube_shape = rss.cube_raw.shape
    shmem_nonlinpoly = multiprocessing.shared_memory.SharedMemory(
        name='nonlinpoly', create=True,
        size=(dummy.itemsize * (poly_order+1) * cube_shape[1] * cube_shape[2]),
    )
    result_nonlinpoly = numpy.ndarray(
        shape=(poly_order+1, cube_shape[1], cube_shape[2]),
        dtype=numpy.float32, buffer=shmem_nonlinpoly.buf)

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
        for x,y in itertools.product(range(2048), range(2048)):

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
                            cube_shape=rss.cube_raw.shape,
                            poly_order=poly_order,
                            ref_level=args.reflevel,
                            saturation_level=args.saturation,
                            workername="Worker_%03d" % (n+1)),
                daemon=True
            )
            jobqueue.put(None)
            p.start()
            worker_processes.append(p)

        # gather results
        logger.info("Gathering results")
        output_cube = numpy.full((poly_order+1,2048,2048), fill_value=numpy.NaN)
        for n in range(n_jobs):
            (x,y,cpu_time,reflevel) = resultqueue.get()
            out_processing_time[y,x] = cpu_time
            out_refslope[y,x] = reflevel
            # output_cube[:,y,x] = polyfit
            # print(polyfit)

        # wait for all work to be done
        logger.info("Working for parallel fitting to complete")
        jobqueue.join()

        # make sure all processes are shut down
        for p in worker_processes:
            p.join()

        if (args.nonlinearity_fn is None):
            out_fn = "nonlinpoly.fits"
        else:
            out_fn = args.nonlinearity_fn
        logger.info("Writing correction coefficients to output FITS (%s)" % (out_fn))
        pyfits.PrimaryHDU(data=result_nonlinpoly).writeto(out_fn, overwrite=True)

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
    shmem_nonlinpoly.close()
    shmem_nonlinpoly.unlink()

    # rss.plot_pixel_curve(818,1033)
    # rss.plot_pixel_curve(1700,555)



if __name__ == "__main__":
    main()