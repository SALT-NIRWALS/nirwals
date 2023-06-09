#!/usr/bin/env python3


import os
import sys
import numpy
import nirwals_reduce

import astropy.io.fits as pyfits

if __name__ == "__main__":
    fn = sys.argv[1]

    print("Loading all files, preparing for plotting")
    rss = rss_reduce.RSS(fn=fn)
    rss.load_all_files()

    amp1s = []
    lefts = []
    mains = []
    rights = []
    for frame in range(rss.image_stack.shape[0]):

        img = rss.image_stack[frame]
        left = numpy.mean(img[:, 1:4], axis=1)
        main = numpy.mean(img[:, 5:-5], axis=1)
        right = numpy.mean(img[:, -4:-1], axis=1)
        amp1 = numpy.mean(img[:,6:62], axis=1)
        # print(left.shape, main.shape)

        lefts.append(left)
        mains.append(main)
        rights.append(right)
        amp1s.append(amp1)
        #break

    lefts = numpy.array(lefts)
    mains = numpy.array(mains)
    rights = numpy.array(rights)
    amp1s = numpy.array(amp1s)
    # print(lefts.shape)

    out_fn = sys.argv[2]
    with open(out_fn, "w") as dumpfile:

        for row in range(lefts.shape[1]):
            print("\n\n\n# row %d" % (row), file=dumpfile)
            numpy.savetxt(dumpfile,
                          numpy.array([lefts[:, row], mains[:, row],
                                       rights[:, row], amp1s[:, row]]).T )
