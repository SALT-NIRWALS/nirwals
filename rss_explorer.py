#!/usr/bin/env python3

import matplotlib
print(matplotlib.get_backend())
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
import threading
import queue

import pyds9
import multiprocessing

import rss_reduce




def ds9_listener(ds9, return_queue):

    while(True):

        # wait for user-interaction from ds9
        try:
            # reply = ds9.get("imexam coordinate image")
            reply = ds9.get("iexam any coordinate image")
            # reply = input("ds9 result")
        except ValueError:
            return_queue.put(None)

            print("Shutting down")
            break

        # interpret the interaction
        try:
            _items = reply.split()
            print(_items)
            command = _items[0]
            ix = int(round(float(_items[1])))-1
            iy = int(round(float(_items[2])))-1
            print(ix,iy)
        except:
            continue

        # forward the answer to the worker_queue for plotting
        return_queue.put((command, ix, iy))


# adopting solution asked in
# https://stackoverflow.com/questions/61397176/
# .. how-to-keep-matplotlib-from-stealing-focus
# and taken from here:
# https://stackoverflow.com/questions/45729092/
# .. make-interactive-matplotlib-window-not-pop-to-front-on-each-update-windows-7/45734500#45734500
def mypause(interval):
    backend = plt.rcParams['backend']
    if backend in matplotlib.rcsetup.interactive_bk:
        figManager = matplotlib._pylab_helpers.Gcf.get_active()
        if figManager is not None:
            canvas = figManager.canvas
            if canvas.figure.stale:
                canvas.draw()
            canvas.start_event_loop(interval)
            return

def rss_plotter(rss, ds9_queue):

    print("Plotter thread running")

    plt.ion()
    fig = plt.figure()
    fig.show()
    plt.show(block=False)

    ax = fig.add_subplot(111)
    plt.pause(0.01)

    while (True):

        print("Waiting for command")
        ds9_command = ds9_queue.get()
        if (ds9_command is None):
            print("Shutting down plotter")
            break

        print("making plot")
        command, ix, iy = ds9_command

        # x = numpy.linspace(0, ix, 100)
        # y = x * iy
        # fig.clf()
        # ax = fig.add_subplot(111)
        # ax.scatter(x,y)
        # # fig.draw()
        # fig.canvas.draw_idle()
        # plt.pause(0.05)
        #
        # time.sleep(5)
        # continue

        ax.cla()
        ax = fig.add_subplot(111)
        fig.suptitle("Pixel position: x=%d // y=%d" % (ix+1, iy+1))

        if (command == "w"):

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
            diff_flux = numpy.pad(numpy.diff(linearized), (1, 0), mode='constant', constant_values=0)
            diff_time = numpy.pad(numpy.diff(rss.read_times), (1, 0), mode='constant', constant_values=0)

            # ax.cla()
            ax.scatter(rss.read_times, diff_flux / diff_time)
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

        # fig.show()
        # plt.show()
        # plt.pause(0.05)
        mypause(0.05)
        # fig.canvas.draw_idle()
        # time.sleep(0.1)

        # ds9.set("cursor "+reply)
    print("Plotter thread shutting down")


if __name__ == "__main__":
    fn = sys.argv[1]

    print("Starting ds9 and establishing connection")
    ds9 = pyds9.DS9() #target='DS9:RSS_Explorer', start=True)

    plt.ion()
    print("Interactive?", plt.isinteractive())

    rss = None
    print("Preparing RSS data cube")
    rss = rss_reduce.RSS(fn=fn, max_number_files=20, use_reference_pixels=True)
    rss.reduce()

    # load image into ds9
    ds9.set_np2arr(rss.weighted_mean)

    ds9_queue = multiprocessing.Queue()

    print("starting ds9 listener thread")
    ds9_thread = multiprocessing.Process(
        target=ds9_listener,
        kwargs=dict(ds9=ds9,
                    return_queue=ds9_queue,)
    )
    ds9_thread.daemon = True
    ds9_thread.start()

    # rss_plotter(rss, ds9_queue)

    print("Starting plotter thread")
    plotter = multiprocessing.Process(
        target=rss_plotter,
        kwargs=dict(rss=rss, ds9_queue=ds9_queue),
    )
    plotter.daemon = True
    plotter.start()


    ds9_thread.join()
    plotter.join()


