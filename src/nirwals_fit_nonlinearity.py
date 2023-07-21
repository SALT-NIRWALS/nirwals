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

            nonlin_results = numpy.polyfit(reads_refpixelcorr[good4fit], (times*slope_reflevel)[good4fit],
                                           deg=poly_order, full=True)
            nonlin_bestfit = nonlin_results[0]
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
    cmdline.add_argument("--saturation", dest="saturation", default=62000,
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

    # rss.reduce(write_dumps=False)
    # rss.write_results()

    rss.load_all_files()
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
            reads_raw = rss.image_stack[:, y,x]
            reads_refpixelcorr = rss.image_stack[:, y,x]
            # print(reads_raw.shape, reads_refpixelcorr.shape)

            jobqueue.put((x, y, reads_raw, reads_refpixelcorr))
            n_jobs += 1

            # if (n_jobs > 100000):
            #     break
            # break
        logger.info("Done with filling queue")
        print("STACK: %d" % (rss.image_stack.shape[0]))

        worker_processes = []
        for n in range(args.n_cores):
            p = multiprocessing.Process(
                target= nonlinfit_worker,
                kwargs=dict(jobqueue=jobqueue,
                            resultqueue=resultqueue,
                            times=rss.raw_read_times[0,:rss.image_stack.shape[0]],
                            poly_order=poly_order,
                            ref_level=10000,
                            saturation_level=55000,
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

                fig = plt.figure()
                ax = fig.add_subplot(111)

                raw_sequence = rss.image_stack[:,y,x]
                raw0 = raw_sequence - numpy.nanmin(raw_sequence)
                read_number = numpy.arange(raw_sequence.shape[0])

                poly = coeffs[:, y,x]

                times = rss.raw_read_times[0,:]
                corrected = numpy.polyval(poly, raw0)
                # print(raw_sequence.shape)
                ax.scatter(times, raw0, s=2, alpha=0.2, c='blue', label='raw')
                ax.plot(times, raw0, 'b-', linewidth=1, c='blue')

                ax.scatter(times, corrected, c='orange', label='linearized', s=1)

                ax.legend(loc='upper left')

                maxy = numpy.min([numpy.max(raw0), 65000])
                maxt = numpy.nanmax(times)
                ax.set_ylim((-0.03*maxy,1.04*maxy))
                ax.set_xlim((-0.03*maxt,1.03*maxt))
                fig.suptitle("Pixel %d , %d" % (x,y))

                pdf.savefig(fig)
                # fig.savefig()
                #

    # rss.plot_pixel_curve(818,1033)
    # rss.plot_pixel_curve(1700,555)
