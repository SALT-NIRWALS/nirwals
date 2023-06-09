#!/usr/bin/env python3


import astropy.io.fits as pyfits
import numpy
import os

import logging
import multiparlog as mplog
import argparse


def read_at_filelist(listfn, check_fn=True):
    filelist = []
    with open(listfn) as lf:
        lines = lf.readlines()
        for line in lines:
            if (line.startswith("#") or len(line.strip()) <= 0):
                continue
            fn = line.split()[0]
            if (check_fn and os.path.isfile(fn)):
                filelist.append(fn)
            else:
                filelist.append(line)

    return filelist


def ramps_by_pixel(filelist, pixel_coord_list):

    #
    # Extract all pixel values, and timings
    #  - Figure out all files for each sequence
    #  - load all input files, apply reference pixel corrections
    #  - extract sequence data for every requested pixel
    #
    exptime_sequences = []
    for fn in filelist:

        bn = ".".join(fn.split(".")[:-2])
        dir,local_bn = os.path.split(bn)
        print(bn)

        # open file, get number of reads/ramps/etc
        ref_hdu = pyfits.open(fn)
        n_reads = ref_hdu[0].header['NGROUPS']
        print(n_reads)

        actexp_list = []
        pixel_fluxes = numpy.full((n_reads, n_pixels), fill_value=numpy.NaN)
        for r in range(n_reads):

            fitsfn = "%s.%d.fits" % (bn, r+1)
            try:
                hdulist = pyfits.open(fitsfn)
            except OSError:
                actexp_list.append(numpy.NaN)
                print("Failed to open file %s " % (fitsfn))
                continue

            # track integration times
            actexp = hdulist[0].header['ACTEXP'] / 1e6
            actexp_list.append(actexp)
            print(fitsfn, actexp)

            # load image, correct reference pixels, and extract pixel fluxes
            img = hdulist[0].data.astype(numpy.float)
            refpixcorr = rss_refpixel_calibrate.reference_pixels_to_background_correction(
                data=img,
                debug=False
            )
            img -= refpixcorr

            for ip, (x,y) in enumerate(pixel_coord_list):
                pixel_fluxes[r,ip] = img[y-1,x-1]

        actexp_list = numpy.array(actexp_list)
        txt_fn = "%s__%s.txt" % (local_bn, args.prefix)
        print(actexp_list.shape, pixel_fluxes.shape)
        numpy.savetxt(txt_fn, numpy.hstack((actexp_list.reshape((-1,1)), pixel_fluxes)))
        print(txt_fn)

        # subtract the average of the three lowest samples as an additional reference
        pixel_fluxes_sorted = numpy.sort(pixel_fluxes, axis=0)
        pixel_fluxes_offset = numpy.mean(pixel_fluxes_sorted[:3, :], axis=0).reshape((1, -1))
        print("offset shape:", pixel_fluxes_offset.shape)
        pixel_fluxes_offset_corrected = pixel_fluxes - pixel_fluxes_offset
        txt_fn = "%s__%s-corr.txt" % (local_bn, args.prefix)
        print(actexp_list.shape, pixel_fluxes_offset_corrected.shape)
        numpy.savetxt(txt_fn, numpy.hstack((actexp_list.reshape((-1, 1)), pixel_fluxes_offset_corrected)))

        # # add some sequence fitting here
        # fit_results = scipy.optimize.least_squares(
        #     fun=_persistency_plus_signal_fit_err_fct,
        #     x0=pinit,
        #     bounds=([-10, 0, 0.2], [numpy.Inf, 65e3, 100]),
        #     kwargs=dict(read_time=read_time, rate=rate, uncert=uncert),
        # )
        # bestfit = fit_results.x
        # bounds_limited_mask = fit_results.active_mask
        # fit_successful = fit_results.success
        # break


def ramp_by_file(filelist, pixellist):
    pass


if __name__ == "__main__":

    mplog.setup_logging(debug_filename="../debug.log",
                        log_filename="../compare_ramps.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("CompareRamps")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--maxfiles", dest="max_number_files", default=None, type=int,
                         help="limit number of files to load for processing")
    cmdline.add_argument("--nonlinearity", dest="nonlinearity_fn", type=str, default=None,
                         help="non-linearity correction coefficients (3-d FITS cube)")
    cmdline.add_argument("--prefix", dest="prefix", type=str, default='sequence',
                         help="non-linearity correction coefficients (3-d FITS cube)")
    # cmdline.add_argument("--persistency", dest="persistency_mode", type=str, default="quick",
    #                      help="persistency mode")
    # cmdline.add_argument("--saturation", dest="saturation", default=62000,
    #                      help="saturation value/file")
    # cmdline.add_argument("healpix32", nargs="+", type=int,
    #                      help="list of input filenames")
    # cmdline.add_argument("--rerun", type=int, default=6,
    #                      help="rerun")
    # cmdline.add_argument("--dumps", dest="write_dumps", default=None,
    #                      help="write intermediate process data [default: NO]")
    # cmdline.add_argument("--debugpngs", dest="write_debug_pngs", default=False, action='store_true',
    #                      help="generate debug plots for all pixels with persistency [default: NO]")
    cmdline.add_argument("--refpixel", dest="use_ref_pixels", default=False, action='store_true',
                         help="use reference pixels [default: NO]")
    # cmdline.add_argument("--flat4salt", dest="write_flat_for_salt", default=False, action='store_true',
    #                      help="write a flat, 1-extension FITS file for SALT")
    cmdline.add_argument("--perframe", dest="combine_pixels_per_frame", default=False, action='store_true',
                         help="combined multiple pixels per file (rathen than multiple files per pixel)")
    cmdline.add_argument("inputs", nargs="+",
                         help="list of input filenames : pixels")
    args = cmdline.parse_args()


    #
    # Read input options: List of input filenames, and pixel coordinates
    #
    filelist = []
    pixel_coord_list = []
    is_pixel = False
    for i in args.inputs:
        if (i == ":"):
            is_pixel = True
            continue
        if (is_pixel):
            print(i)
            if (i.startswith("@") and os.path.isfile(i[1:])):
                add_list = read_at_filelist(i[1:], check_fn=False)
                print(add_list)
                new_xy = [[int(float(x)) for x in l.split()[:2]] for l in add_list]
                pixel_coord_list.extend(new_xy)
            else:
                xy = [int(float(x)) for x in i.split(",")[:2]]
                pixel_coord_list.append(xy)
        else:
            if (i.startswith("@") and os.path.isfile(i[1:])):
                filelist.extend(read_at_filelist(i[1:]))
            elif (os.path.isfile(i)):
                filelist.append(i)

    print("Files:", "\n -- ".join(['']+filelist))
    print("Pixels:", pixel_coord_list)
    n_pixels = len(pixel_coord_list)

    ramps_by_pixel(filelist, pixel_coord_list)
