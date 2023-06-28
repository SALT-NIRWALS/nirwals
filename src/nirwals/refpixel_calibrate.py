#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import numpy

import matplotlib.pyplot as plt
import logging
import multiparlog as mplog
import argparse

dummycounter = 0

def reference_pixels_to_background_correction(data, edge=1, verbose=False, make_plots=False, debug=False, yslope=True, even_odd=False):

    global dummycounter
    dummycounter += 1

    # # first, combine left & right to subtract row-wise overscan level
    # _left = numpy.mean(data[:, edge:4], axis=1).reshape((-1,1))
    # _right = numpy.mean(data[:, -4:-edge], axis=1).reshape((-1,1))
    # row_wise = numpy.mean([_left, _right], axis=0)
    # if (debug):
    #     print(row_wise.shape)

    # plt.scatter(numpy.arange(row_wise.shape[0]), row_wise, s=1)
    # plt.show()

    # data_rowsub = data - row_wise
    # if (debug):
    #     pyfits.PrimaryHDU(data=data_rowsub).writeto("del__rowsub_%d.fits" % (dummycounter), overwrite=True)

    # now figure out the column situation
    top = data[edge:4, :]
    bottom = data[-4:-edge, :]
    ref_columns = numpy.vstack([top, bottom])
    if (debug):
        print("ref-columns shape:", ref_columns.shape)
    n_rows = 8 - 2*edge
    n_amps = 32
    amp_size = ref_columns.shape[1] // n_amps

    # create some fake signal to check if the folding works as expected
    iy,ix = numpy.indices(ref_columns.shape)
    # ref_columns = numpy.cos(ix * numpy.pi / 32) + 0.5*numpy.cos(ix * numpy.pi/64)

    if (yslope):
        # If we are fitting vertical slopes, this is done on a column-by-column basis
        mean_top = numpy.mean(top, axis=0)
        mean_bottom = numpy.mean(bottom, axis=0)
        slopes = (mean_bottom - mean_top) / (2048-edge-4)  # that's the number of pixels between top & bottom

        full_iy,_ = numpy.indices(data.shape, dtype=float)
        full_2d_correction = full_iy * slopes + mean_top

        if (debug):
            pyfits.PrimaryHDU(data=full_2d_correction).writeto("refpixel__yslope.fits", overwrite=True)

        return full_2d_correction

    else:

        # take out the average intensity in each amplifier-block
        amp_blocks = ref_columns.T.reshape((n_amps, -1)) #.reshape((-1, n_amps))
        avg_amp_background = numpy.median(amp_blocks, axis=1)
        amp_background = numpy.repeat(avg_amp_background.reshape((-1,1)), amp_size).reshape((1,-1))

        # prepare the first-order correction and apply to the reference pixels themselves
        ref_background = numpy.median(ref_columns, axis=0)
        amp_background_2d = numpy.ones_like(ref_columns) * amp_background

        if (debug):
            print("amp_blocks shape:", amp_blocks.shape)
            print("Avg amp levels:", avg_amp_background.shape)
            print("amp background:", amp_background.shape)

            print("REF BG:", ref_background.shape)
            numpy.savetxt("ref_cols.txt", ref_background)
            numpy.savetxt("amp_background.txt", amp_background.T)
            numpy.savetxt("ref_columns.txt", ref_columns.T)

    if (make_plots):
        fig = plt.figure()
        fig.suptitle("ref-pixel average intensity -- data vs median")
        ax = fig.add_subplot(111)
        ax.scatter(numpy.arange(ref_columns.shape[1]), ref_columns[0, :], s=0.2)
        ax.scatter(numpy.arange(ref_columns.shape[1]), ref_columns[2, :], s=0.2)
        ax.scatter(numpy.arange(ref_columns.shape[1]), ref_columns[5, :], s=0.2)
        ax.scatter((numpy.arange(n_amps)+0.5)*amp_size, avg_amp_background)
        fig.show()


    normalized_ref_columns = ref_columns - amp_background

    if (not even_odd):
        # we only correct for the bias level

        print("normalized ref columns shape:", normalized_ref_columns.shape)
        numpy.savetxt("norm_ref_columns.txt", normalized_ref_columns.T)
        print("amp_background_2d shape", amp_background_2d.shape)

        full_2d_correction = numpy.ones_like(data) * amp_background
        print("full 2d corr shape", full_2d_correction.shape)
        return full_2d_correction



    # reshape the ref columns to align pixels read out in parallel
    if (debug):
        print("every amp has %d pixels" % (amp_size))
    ref_cols_2amps = normalized_ref_columns.reshape(-1, 2*amp_size)

    # now flip one half to account for different read-directions
    ref_cols_2amps_flipped = numpy.array(ref_cols_2amps)
    ref_cols_2amps_flipped[:, amp_size:] = numpy.flip(ref_cols_2amps[:, amp_size:], axis=1)
    if (debug):
        print("2amps:", ref_cols_2amps.shape)

    # align all pixels read out simultaneous
    ref_cols_1amp = ref_cols_2amps_flipped.reshape(-1, amp_size)
    if (debug):
        print("1amp:", ref_cols_1amp.shape)

    # calculate correction
    ref_cols_combined = numpy.mean(ref_cols_1amp, axis=0)
    if (debug):
        print("combined signal:", ref_cols_combined.shape)

    # with this correction, reconstruct the full column-wise correction
    ref_cols_correction_2amp = numpy.hstack([ref_cols_combined, numpy.flip(ref_cols_combined)]).reshape((1,-1))
    ref_cols_correction_full = numpy.repeat(ref_cols_correction_2amp, 16, axis=0).reshape((1,-1))
    if (debug):
        print("full correction:", ref_cols_correction_full.shape)

    ref_cols_correction_full_2d = numpy.repeat(ref_cols_correction_full, 10, axis=0)
    if (debug):
        print("2d correction:", ref_cols_correction_full_2d.shape)

    total_column_correction = ref_cols_correction_full + amp_background
    if (debug):
        print("total col correction shape:", total_column_correction.shape)
        numpy.savetxt("total_column_corr.txt", total_column_correction[0,:])


    # apply the column-wise correction to the full frame
    image = data - total_column_correction

    if (debug):
        pyfits.PrimaryHDU(data=image).writeto("del__totalcorrected.fits", overwrite=True)

    full_2d_correction = total_column_correction
    if (debug):
        print("full 2d:", full_2d_correction.shape)

    # combine
    if (debug):
        pyfits.HDUList([
            pyfits.PrimaryHDU(),
            pyfits.ImageHDU(data=ref_columns, name="REF_COLS"),
            pyfits.ImageHDU(data=amp_background_2d, name="REF_COLS_background"),
            pyfits.ImageHDU(data=normalized_ref_columns, name="NORM_REF_COLS"),
            pyfits.ImageHDU(data=ref_cols_2amps, name="DOUBLE_AMP"),
            pyfits.ImageHDU(data=ref_cols_2amps_flipped, name="DOUBLE_AMP_FLIPPED"),
            pyfits.ImageHDU(data=ref_cols_1amp, name="SINGLE_AMP"),
            pyfits.ImageHDU(data=ref_cols_correction_2amp, name="CORR_2AMP"),
            pyfits.ImageHDU(data=ref_cols_correction_full_2d, name="CORR_2D"),
            pyfits.ImageHDU(data=full_2d_correction, name="TOTAL_CORR"),
        ]).writeto("del__refcols.fits", overwrite=True)


    return full_2d_correction



if __name__ == "__main__":




    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("NIRWALS-RefPixel")

    cmdline = argparse.ArgumentParser()
    cmdline.add_argument("--output", dest="output_fn", type=str, default="reduced",
                         help="addition to output filename")
    cmdline.add_argument("--evenodd", dest="even_odd", default=False, action='store_true',
                         help="apply even/odd pixel correction")
    cmdline.add_argument("--yslope", dest="yslope", default=False, action='store_true',
                         help="apply vertical slope correction")
    cmdline.add_argument("files", nargs="+",
                         help="list of input filenames")
    args = cmdline.parse_args()

    fn = args.files[0] #sys.argv[1]
    output_fn = args.output_fn #sys.argv[2]

    hdulist = pyfits.open(fn)
    data = hdulist[0].data

    full_2d_correction = reference_pixels_to_background_correction(
        data, debug=True, make_plots=True,
        even_odd=args.even_odd, yslope=args.yslope)

    data = data - full_2d_correction
    hdulist[0].data = data

    hdulist.writeto(output_fn, overwrite=True)

