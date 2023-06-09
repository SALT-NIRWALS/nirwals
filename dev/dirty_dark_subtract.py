#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits as pyfits


if __name__ == "__main__":

    sci_fn = sys.argv[1]
    dark_fn = sys.argv[2]

    out_fn = sys.argv[3]

    sci_hdu = pyfits.open(sci_fn)
    sci = sci_hdu['PERS.SIGNAL'].data

    dark_hdu = pyfits.open(dark_fn)
    dark = dark_hdu['PERS.SIGNAL'].data

    darksub = sci - dark

    out_hdu = pyfits.PrimaryHDU(data=darksub)
    out_hdu.header['SCI'] = sci_fn
    out_hdu.header['DARK'] = dark_fn
    out_hdu.writeto(out_fn, overwrite=True)
    