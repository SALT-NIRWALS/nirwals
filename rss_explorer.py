#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits
import matplotlib.pyplot as plt
import time

import pyds9

import rss_reduce

plt.ion()


if __name__ == "__main__":

    fn = sys.argv[1]

    print("Starting ds9 and establishing connection")
    ds9 = pyds9.DS9() #target='DS9:RSS_Explorer', start=True)

    print("Preparing RSS data cube")
    rss = rss_reduce.RSS(fn=fn, max_number_files=20, use_reference_pixels=True)
    rss.reduce()

    ds9.set_np2arr(rss.weighted_mean)

    fig = plt.figure()
    fig.show()
    ax = fig.add_subplot(111)
    while (True):

        try:
            reply = ds9.get("imexam coordinate image")
        except ValueError:
            print("Shutting down")
            break

        # reply = ds9.get("crosshair image")
        # print(reply)

        try:
            _items = reply.split()
            print(_items)
            ix = int(round(float(_items[0])))-1
            iy = int(round(float(_items[1])))-1
            print(ix,iy)
        except:
            continue

        ax.cla()
        ax.scatter(rss.read_times, rss.image_stack[:, iy, ix])
        fig.show()

        # ds9.set("cursor "+reply)


