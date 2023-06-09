#!/usr/bin/env python3

import os
import sys
import numpy
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt

import nirwals_reduce

if __name__ == "__main__":

    fn = sys.argv[1]
    rss = rss_reduce.RSS(fn)
    rss.load_all_files()

    # rss.read_nonlinearity_corrections("nonlin_inverse.fits")
    rss.read_nonlinearity_corrections("nonlin3d.fits")

    coord_list_fn = sys.argv[2]
    lines = []
    with open(coord_list_fn, "r") as f:
        lines = f.readlines()
    for l in lines:
        if (l.startswith("#")):
            continue
        items = l.split()
        print(items)
        if (len(items) <= 2):
            continue
        x = int(float(items[0]))
        y = int(float(items[1]))
        #x = int(numpy.round(float(items[0])),0) - 1
        #y = int(numpy.round(float(items[1])),0) - 1

        stack = rss.image_stack[:, y,x]

        stack_zeroed = stack - stack[0]

        factors = rss.nonlinearity_cube[:, y,x]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        _i = numpy.arange(stack.shape[0])

        slope = stack_zeroed[10] / 10
        linear = slope * _i

        corrected = factors[0] * stack_zeroed + \
                    factors[1] * numpy.power(stack_zeroed, 2) + \
                    factors[2] * numpy.power(stack_zeroed, 3)

        ax.scatter(_i, stack_zeroed)
        ax.scatter(_i, corrected)
        ax.plot(_i, corrected)
        ax.plot(_i, linear)
        fig.savefig("%s___nonlin_%04d_%04d.png" % (rss.filebase, x,y))
        plt.close(fig)

