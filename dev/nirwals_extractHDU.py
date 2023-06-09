#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits


if __name__ == "__main__":

    in_fn = sys.argv[1]
    extensions = sys.argv[2:-1]
    out_fn = sys.argv[-1]

    hdulist = pyfits.open(in_fn)
    out_list = [hdulist[0]]
    for ext in extensions:
        # try:
        out_list.append(hdulist[ext])
        # except:
        #     print("Unable to find extension %s" % (ext))


    out_hdulist = pyfits.HDUList(out_list)
    out_hdulist.writeto(out_fn, overwrite=True)