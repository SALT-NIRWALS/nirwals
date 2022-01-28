#!/usr/bin/env python3


import os
import sys
import numpy

import astropy.io.fits as pyfits

if __name__ == "__main__":

    fn = sys.argv[1]

    hdulist = pyfits.open(fn)
    data = hdulist[0].data.astype(numpy.float32)

    # extract top and bottom N lines
    n_lines = 4
    n_amplifiers = 32

    top = data[:n_lines, :]
    bottom = data[-n_lines:, :]
    top_and_bottom = numpy.vstack([top, bottom])
    print(top_and_bottom.shape)

    # combine all columns
    flattened = numpy.mean(top_and_bottom, axis=0)
    print(flattened.shape)

    # now reshape into even/odd pairs
    even_odd = flattened.reshape((n_amplifiers*2, -1))
    print(even_odd.shape)

    # combine all even/odd pairs
    even_odd_levels = numpy.mean(even_odd, axis=1)
    print(even_odd_levels)
    print(even_odd_levels.shape)

    # duplicate back to full-res
    # overscan_levels = numpy.repeat(even_odd_levels.reshape((-1,2)), even_odd.shape[1], axis=1).reshape((1,-1))
    overscan_levels = even_odd_levels.reshape((-1,2)).repeat(even_odd.shape[1], axis=0).reshape((1,-1))
    print(overscan_levels[:5])
    print("overscan_levels:", overscan_levels.shape)

    subtracted = data - overscan_levels

    all_overscan = numpy.zeros_like(data) + overscan_levels

    # hdulist[0].data = data

    combined = data * 1.0
    combined[4:10, :] = all_overscan[4:10, :]

    out_fn = sys.argv[2] #"even_odd.fits"
    print("Writing results to %s" % (out_fn))
    out_hdu = pyfits.HDUList([
        pyfits.PrimaryHDU(header=hdulist[0].header),
        pyfits.ImageHDU(data=data, name="RAW"),
        pyfits.ImageHDU(data=subtracted, name="SUBTRACTED"),
        pyfits.ImageHDU(data=all_overscan, name="OVERSCAN"),
        pyfits.ImageHDU(data=combined, name="COMBINED"),
    ])
    out_hdu.writeto(out_fn, overwrite=True)





