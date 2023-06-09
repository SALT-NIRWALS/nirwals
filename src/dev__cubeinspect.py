#!/usr/bin/env python3


import os
import sys
import numpy

import nirwals_reduce
import pyds9 #as pyds9
import matplotlib.pyplot as plt
import time

if __name__ == "__main__":

    fn = sys.argv[1]

    display_fn = sys.argv[2]

    ds9 = pyds9.DS9(start=True) #, wait=10, verify=True)
    print(ds9)

    print("Loading display frame into ds9")
    ds9.set("file %s" % (display_fn))

    print("Loading all files, preparing for plotting")
    rss = rss_reduce.RSS(fn=fn)
    rss.load_all_files()

    while (True):
        pos = ds9.get("imexam coordinate image")
        # print(pos, type(pos))

        [x,y] = [int(numpy.round(float(i),0)) for i in pos.split()]
        print(pos, type(pos), "-->", x, y)
        # time.sleep(1)

        # fig = plt.figure()

        rss.plot_pixel_curve(x,y, differential=True, diff_vs_cum=False, cumulative=True,
                             show_plot=True)
