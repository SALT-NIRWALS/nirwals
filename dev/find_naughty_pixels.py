#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import numpy


if __name__ == "__main__":

    fn = sys.argv[1]
    hdulist = pyfits.open(fn)

    median = hdulist['MEDIAN'].data
    sigma = hdulist['SIGMA'].data

    gain = 1
    readnoise = 20

    noise_expected = numpy.sqrt( median * gain + 2 * readnoise ** 2) / gain

    noise_relative = sigma / noise_expected

    out_hdulist = pyfits.HDUList([
        pyfits.PrimaryHDU(),
        hdulist['MEDIAN'],
        hdulist['SIGMA'],
        pyfits.ImageHDU(data=noise_expected, name="NOISE_EXPECTED"),
        pyfits.ImageHDU(data=noise_relative, name='NOISE_RELATIVE'),
    ])

    out_fn = fn[:-5]+"__complete.fits"
    print("writing %s" % (out_fn))
    out_hdulist.writeto(out_fn, overwrite=True)

    bad_pixels = noise_relative > 5
    bpm_fn = fn[:-5]+"__badpixels.fits"
    pyfits.PrimaryHDU(data=bad_pixels.astype(numpy.int8)).writeto(bpm_fn, overwrite=True)
