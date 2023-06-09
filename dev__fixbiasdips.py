#!/usr/bin/env python3


import os
import sys
import numpy
import nirwals_reduce

import astropy.io.fits as pyfits

if __name__ == "__main__":

    fn = sys.argv[1]

    print("Loading all files, preparing for plotting")
    #rss = rss_reduce.RSS(fn=fn)
    #rss.load_all_files()

    hdulist = pyfits.open(fn)
    data = hdulist[0].data

    n_amps = 32
    amp_width = data.shape[1]//n_amps
    out_data = numpy.array(data, dtype=numpy.float)
    out_data[:, 0:4] = numpy.NaN
    out_data[:, -4:] = numpy.NaN
    for n in range(1, n_amps, 2):

        x1 = n*amp_width
        x2 = x1 + amp_width

        strip = out_data[:, x1:x2]
        flipped = strip[:, ::-1]

        out_data[:, x1:x2] = flipped

    amp_stacked = numpy.nanmean(out_data.reshape((-1, n_amps, amp_width)), axis=1)
    print(amp_stacked.shape)


    out_fn = sys.argv[2]
    out_hdulist = pyfits.HDUList([
        pyfits.PrimaryHDU(header=hdulist[0].header),
        pyfits.ImageHDU(data=out_data, name="amp_flipped"),
        pyfits.ImageHDU(data=amp_stacked, name="amp_stacked"),
        ])

    print("saving output to %s" % (out_fn))
    out_hdulist.writeto(out_fn, overwrite=True)

    print("all done!")

