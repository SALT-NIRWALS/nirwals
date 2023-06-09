#!/usr/bin/env python3

import numpy
import sys
import os
import astropy.io.fits as pyfits

subtract_first = True

if __name__ == "__main__":

    out_fn = sys.argv[1]

    filelist = sys.argv[2:]

    imglist = []
    for fn in filelist:

        if (not os.path.isfile(fn)):
            continue
        try:
            hdulist = pyfits.open(fn)
        except:
            print("Trouble reading from %s" % (fn))
            continue

        # load data, apply reference pixel correction
        print("Working on %s" % (fn))
        raw = hdulist[0].data
        refcorr = rss_refpixel_calibrate.reference_pixels_to_background_correction(
            data=raw, verbose=False, make_plots=False, debug=False
        )
        corr = raw - refcorr

        if (subtract_first):
            items = fn.split(".")
            first_fn = ".".join(items[:-2])+".1.fits"
            print(first_fn)

            first_hdu = pyfits.open(first_fn)
            first_raw = first_hdu[0].data
            first_ref = rss_refpixel_calibrate.reference_pixels_to_background_correction(
                data=first_raw, verbose=False, make_plots=False, debug=False
            )
            first_corr = first_raw - first_ref

            corr -= first_corr

        imglist.append(corr)

    # sys.exit(0)

    # combine individual reads into datacube
    print("Forming datacube")
    imgcube = numpy.array(imglist)

    # calculate sigma deviations
    for iter in range(3):
        print("Computing statistics, iteration %d" % (iter+1))

        _stats = numpy.nanpercentile(imgcube, [16,50,84], axis=0)
        print(_stats[1].shape)

        _median = _stats[1]
        _sigma = 0.5 * (_stats[2] - _stats[0])

        bad = (imgcube > (_median + 3*_sigma)) | (imgcube < (_median - 3*_sigma))
        imgcube[bad] = numpy.NaN



        outhdu = pyfits.HDUList([
            pyfits.PrimaryHDU(),
            pyfits.ImageHDU(data=_median, name="MEDIAN"),
            pyfits.ImageHDU(data=_sigma, name="SIGMA")
        ])
        _fn = "%s_%d.fits" % (out_fn, iter+1)
        print("Saving results to %s" % (_fn))
        outhdu.writeto(_fn, overwrite=True)



