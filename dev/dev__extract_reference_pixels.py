#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import numpy


if __name__ == "__main__":

    fn = sys.argv[1]

    hdulist = pyfits.open(fn)
    img = hdulist[0].data

    columns = numpy.zeros((8, img.shape[1]))
    columns[0:4, :] = img[0:4, :]
    columns[-4:, :] = img[-4:, :]

    rows = numpy.zeros((img.shape[0], 8))
    rows[:, 0:4] = img[:, :4]
    rows[:, -4:] = img[:, -4:]

    outhdu = pyfits.HDUList([
        pyfits.PrimaryHDU(),
        pyfits.ImageHDU(data=rows),
        pyfits.ImageHDU(data=columns),
    ])

    outhdu.writeto("reference_only.fits")