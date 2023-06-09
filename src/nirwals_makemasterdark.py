#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits as pyfits

if __name__ == "__main__":

    out_fn = sys.argv[1]
    fn_list = sys.argv[2:]

    dark_stack = []
    for fn in fn_list:

        hdulist = pyfits.open(fn)
        try:
            pers_signal = hdulist['PERS.SIGNAL'].data
            print("Read %s" % (fn))
        except:
            print("Unable to find PERS.SIGNAL in %s" % (fn))
            continue

        # mask out bad pixels
        bad_pixels = (pers_signal < -9)
        pers_signal[bad_pixels] = numpy.NaN

        dark_stack.append(pers_signal)

    dark_stack = numpy.array(dark_stack)

    dark_signal = numpy.nanmean(dark_stack, axis=0)
    dark_rms = numpy.sqrt(numpy.nanvar(dark_stack, axis=0))
    dark_s2n = dark_signal / dark_rms

    dark_masked = dark_signal.copy()
    dark_masked[dark_s2n < 2] = 0.

    print(dark_signal.shape)

    out_hdulist = pyfits.HDUList([
        pyfits.PrimaryHDU(),
        pyfits.ImageHDU(data=dark_signal, name='DARK'),
        pyfits.ImageHDU(data=dark_rms, name='DARK.RMS'),
        pyfits.ImageHDU(data=dark_s2n, name='DARK.S2N'),
        pyfits.ImageHDU(data=dark_masked, name='DARK.MASKED'),
    ])
    print("Writing master-dark to %s" % (out_fn))
    out_hdulist.writeto(out_fn, overwrite=True)

    print("all done!")
