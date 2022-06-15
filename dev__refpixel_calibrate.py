#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import numpy


if __name__ == "__main__":

    fn = sys.argv[1]
    hdulist = pyfits.open(fn)
    data = hdulist[0].data

    edge = 1

    iy,ix = numpy.indices(data.shape, dtype=numpy.float)
    _top = numpy.mean(data[edge:4, :], axis=0).reshape((1,-1))
    _bottom = numpy.mean(data[-4:-edge, :], axis=0).reshape((1,-1))
    print(_top.shape)

    fy = (iy-4)/2040.
    fx = (ix-4)/2040.
    y_pattern = _bottom + (fy * (_top-_bottom))

    y_sub = data - y_pattern
    _left = numpy.mean(y_sub[:, edge:4], axis=1).reshape((-1,1))
    _right = numpy.mean(y_sub[:, -4:-edge], axis=1).reshape((-1,1))
    print(_left.shape)
    x_pattern = _left + fx * (_right-_left)

    pyfits.PrimaryHDU(data=y_pattern).writeto("ypattern.fits", overwrite=True)
    pyfits.PrimaryHDU(data=x_pattern).writeto("xpattern.fits", overwrite=True)

    xy_pattern = y_pattern #+ x_pattern
    pyfits.PrimaryHDU(data=xy_pattern).writeto("xypattern.fits", overwrite=True)

    pattern_sub = data - (y_pattern + x_pattern)
    pyfits.PrimaryHDU(data=pattern_sub).writeto("pattern_sub.fits", overwrite=True)


    zero_read = pyfits.open("210513_LED_1550nm_18dB_dewar_120K_encl_22C.350.1.1.fits")[0].data.astype(numpy.float)
    zerosub = (data - zero_read)
    pyfits.PrimaryHDU(data=zerosub).writeto("zero_sub.fits", overwrite=True)

