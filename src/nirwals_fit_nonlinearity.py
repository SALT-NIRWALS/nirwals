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



def nonlinfit_worker(jobqueue, resultqueue, times, poly_order=3, ref_level=10000, saturation_level=55000, workername="NonLinFitWorker"):

    logger = logging.getLogger(workername)
    logger.info("Starting worker %s" % (workername))

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

        x, y, reads_refpixelcorr, reads_raw = job
        # logger.debug("x=%d, y=%d: read:%s raw:%s times:%s" % (x,y,str(reads_refpixelcorr.shape), str(reads_raw.shape), str(times.shape)))
        # print(times)
        # print(reads_refpixelcorr)
        # print(reads_raw)

        # subtract off any residual offsets (read for t=0 should be 0)
        reads_offset = numpy.nanmin(reads_refpixelcorr)
        reads_refpixelcorr -= reads_offset

        # numpy.savetxt("dummydump_%04d_%04d.deleteme" % (x,y), numpy.array([
        #     times,reads_refpixelcorr,reads_raw
        # ]).T)

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

        resultqueue.put((x,y,nonlin_bestfit,t2-t1))
        jobqueue.task_done()

    logger.info("Shutting down worker %s" % (workername))



if __name__ == "__main__":

    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
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
    rss = NIRWALS(fn, saturation=saturation_fn,
                  max_number_files=args.max_number_files,
                  use_reference_pixels=args.ref_pixel_mode,)

    rss.load_all_files()

    rss.apply_reference_pixel_corrections()
    # rss.reduce(write_dumps=False)
    # rss.write_results()

    # rss.subtract_first_read()

    if (not args.verify):
        # rss.fit_nonlinearity(ref_frame_id=4, make_plot=False)

        jobqueue = multiprocessing.JoinableQueue()
        resultqueue = multiprocessing.Queue()
        times = rss.raw_read_times

        poly_order = 5
        logger.info("Starting to fill queue")
        n_jobs = 0
        for x,y in itertools.product(range(2048), range(2048)):

            # while we use cube_linearized, we did not actually apply any nonlinearity corrections
            reads_raw = rss.cube_raw[:, y,x]
            reads_refpixelcorr = rss.cube_linearized[:, y,x]
            # print(reads_raw.shape, reads_refpixelcorr.shape)

            jobqueue.put((x, y, reads_raw, reads_refpixelcorr))
            n_jobs += 1

            # if (n_jobs > 100000):
            #     break
            # break
        logger.info("Done with filling queue")
        print("STACK: %d" % (rss.cube_raw.shape[0]))

        worker_processes = []
        for n in range(args.n_cores):
            p = multiprocessing.Process(
                target= nonlinfit_worker,
                kwargs=dict(jobqueue=jobqueue,
                            resultqueue=resultqueue,
                            times=rss.raw_read_times[:rss.cube_raw.shape[0]],
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
            (x,y,polyfit,cpu_time) = resultqueue.get()
            output_cube[:,y,x] = polyfit
            # print(polyfit)

        # wait for all work to be done
        logger.info("Working for parallel fitting to complete")
        jobqueue.join()

        # make sure all processes are shut down
        for p in worker_processes:
            p.join()

        out_fn = "nonlinpoly.fits"
        logger.info("Writing correction coefficients to output FITS (%s)" % (out_fn))
        pyfits.PrimaryHDU(data=output_cube).writeto(out_fn, overwrite=True)

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

    # rss.plot_pixel_curve(818,1033)
    # rss.plot_pixel_curve(1700,555)
