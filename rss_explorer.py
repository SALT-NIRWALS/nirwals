#!/usr/bin/env python3

import matplotlib
# matplotlib.use('QT5Agg')
# matplotlib.use('GTK3Agg')
# matplotlib.use('WebAgg')
# matplotlib.use('TkAgg')
# matplotlib.use('WebAgg')


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
    rss = rss_reduce.RSS(fn=fn, max_number_files=50, use_reference_pixels=True)
    rss.reduce()

    ds9.set_np2arr(rss.weighted_mean)

    plt.ion()
    fig = plt.figure()
    fig.show()
    # plt.show()
    # ax = fig.add_subplot(111)
    while (True):

        try:
            # reply = ds9.get("imexam coordinate image")
            reply = ds9.get("iexam any coordinate image")
            # reply = input("ds9 result")
        except ValueError:
            print("Shutting down")
            break

        # reply = ds9.get("crosshair image")
        # print(reply)

        try:
            _items = reply.split()
            print(_items)
            command = _items[0]
            ix = int(round(float(_items[1])))-1
            iy = int(round(float(_items[2])))-1
            print(ix,iy)
        except:
            continue

        # ax.cla()
        ax = fig.add_subplot(111)
        if (command == "w"):
            # ax.cla()

            img_flux = rss.image_stack[:, iy, ix]
            min_flux = numpy.min(img_flux)
            fit_flux = rss.read_times * rss.weighted_mean[iy, ix]

            ax.scatter(rss.read_times, img_flux)
            ax.set_xlabel("Integration time [seconds]")
            ax.set_ylabel("counts")

        elif (command == "s"):
            # ax.cla()

            lin_flux = rss.linearized_cube[:, iy, ix]
            print(lin_flux)
            min_count = numpy.nanmin(lin_flux[1:])
            print(min_count)
            fit_line = rss.read_times * rss.weighted_mean[iy, ix] + min_count

            ax.scatter(rss.read_times, lin_flux)
            ax.plot(rss.read_times, fit_line, "b-")
            ax.set_xlabel("Integration time [seconds]")
            ax.set_ylabel("counts")

        elif (command == "r"):

            linearized = rss.linearized_cube[:, iy, ix]
            diff_flux = numpy.pad(numpy.diff(linearized), (1,0), mode='constant', constant_values=0)
            diff_time = numpy.pad(numpy.diff(rss.read_times), (1,0), mode='constant', constant_values=0)

            # ax.cla()
            ax.scatter(rss.read_times, diff_flux/diff_time)
            ax.axhline(rss.weighted_mean[iy, ix], linestyle='-', color='blue')
            ax.axhline(0., linestyle='--', color="black")
            ax.set_xlabel("Integration time")
            ax.set_ylabel("Differential flux increase between reads [counts/second]")
            ax.set_title("title here")

        elif (command == 'h'):
            print("Stay tuned, help is on the way")

        else:
            print("command (%s) not understood -- press -h- for help")
            continue

        fig.show()
        plt.show()
        # fig.canvas.draw_idle()
        # time.sleep(0.001)

        # ds9.set("cursor "+reply)


