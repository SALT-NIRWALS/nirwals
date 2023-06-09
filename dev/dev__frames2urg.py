#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits as pyfits
import glob


if __name__ == "__main__":

    out_basename = sys.argv[1]

    in_filelist = sys.argv[2:]

    # first, check all files for exposure times



