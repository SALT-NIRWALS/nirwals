#!/usr/bin/env python3

import matplotlib
print(matplotlib.get_backend())
# matplotlib.use('QT5Agg')
# matplotlib.use('GTK3Agg')
# matplotlib.use('WebAgg')
matplotlib.use('TkAgg')
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

import multiparlog as mplog
import logging


import nirwals_reduce




def ds9_listener(ds9, return_queue):

    while(True):

        # wait for user-interaction from ds9
        try:
            # reply = ds9.get("imexam coordinate image")
            reply = ds9.get("iexam any coordinate image")
            # reply = input("ds9 result")
        except TypeError as te:
            print(te)
            continue

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

def rss_plotter(rss, ds9_queue, save_plots=False):

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

        # ax.cla()
        fig.clf()
        ax = fig.add_subplot(111)
        fig.suptitle("Pixel position: x=%d // y=%d" % (ix+1, iy+1))

        # print(rss.read_times)

        pixel_noise = numpy.sqrt(rss.linearized_cube[:, iy, ix])

        if (command == "w"):

            img_flux = rss.image_stack[:, iy, ix]
            min_flux = numpy.min(img_flux)
            fit_flux = rss.read_times * rss.weighted_mean[iy, ix]

            ax.scatter(rss.read_times, img_flux)
            ax.set_xlabel("Integration time [seconds]")
            ax.set_ylabel("counts")
            ax.set_title("Raw read counts [no corrections]")

            if (save_plots):
                plot_filename = "%s___%04dx%04d___rawreads.png" % (rss.filebase, ix+1, iy+1)
                fig.savefig(plot_filename, tight_layout=True, dpi=200)

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
            ax.set_title("Read counts [ref pixel corrected & linearized]")

            if (save_plots):
                plot_filename = "%s___%04dx%04d___netreads.png" % (rss.filebase, ix+1, iy+1)
                fig.savefig(plot_filename, tight_layout=True, dpi=200)

        elif (command == "r" or command=="t"):

            linearized = rss.linearized_cube[:, iy, ix]
            # diff_flux = numpy.pad(numpy.diff(linearized), (1, 0), mode='constant', constant_values=0)
            diff_flux = rss.differential_cube[:, iy, ix]
            diff_time = numpy.pad(numpy.diff(rss.read_times), (1, 0), mode='constant', constant_values=0)

            max_flux = numpy.nanmax(diff_flux + 2*pixel_noise)
            min_flux = numpy.nanmin(diff_flux - 2*pixel_noise)

            # ax.cla()
            ax.scatter(rss.read_times, diff_flux) # / diff_time)
            ax.axhline(rss.weighted_mean[iy, ix], linestyle='-', color='blue')

            if (command == "r"):
                max_flux = numpy.nanmax(diff_flux + 2 * pixel_noise)
                min_flux = numpy.nanmin(diff_flux - 2 * pixel_noise)
                ax.set_ylim(min_flux, max_flux)
            else:
                n_start = diff_flux.shape[0] * 2 // 3
                max_flux = numpy.nanmax(diff_flux[n_start:])
                min_flux = numpy.nanmin(diff_flux[n_start:])
                buf = (max_flux - min_flux) / 10.
                ax.set_ylim(min_flux-buf, max_flux+buf)
            try:
                # also plot the persistency fit, if available
                print("Checking for persistency")
                persistency_fit_data = rss.persistency_fit_global[:, iy, ix]
                print(ix, iy, persistency_fit_data)
                pf = rss_reduce._persistency_plus_signal_fit_fct(persistency_fit_data[:3], rss.read_times)
                ax.plot(rss.read_times, pf)
                ax.set_title("fit flux: %.3f cts/s [%.2f] (persistency: %.2g / tau=%.2f s)" % (
                    persistency_fit_data[0],  rss.weighted_mean[iy, ix],
                    persistency_fit_data[1], persistency_fit_data[2]
                ))
            except:
                ax.set_title("avg flux: %.2f cts/s" % (rss.weighted_mean[iy, ix]))
                pass
            ax.axhline(0., linestyle='--', color="black")
            ax.set_xlabel("Integration time")
            ax.set_ylabel("Differential flux increase between reads [counts/second]")

            if (save_plots):
                plot_filename = "%s___%04dx%04d___rate.png" % (rss.filebase, ix+1, iy+1)
                fig.savefig(plot_filename, tight_layout=True, dpi=200)

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

    mplog.setup_logging(debug_filename="../../debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)


    fn = sys.argv[1]
    try:
        persistency_fn = sys.argv[2]
    except:
        persistency_fn = None

    display_filename = None
    try:
        display_filename = sys.argv[3]
    except:
        pass

    print("Starting ds9 and establishing connection")
    ds9 = pyds9.DS9() #target='DS9:RSS_Explorer', start=True)

    plt.ion()
    print("Interactive?", plt.isinteractive())

    rss = None
    print("Preparing RSS data cube")
    rss = rss_reduce.RSS(fn=fn, max_number_files=50, use_reference_pixels=True)
    rss.reduce()
    print("Using persistency from %s" % (persistency_fn))
    rss.load_precalculated_results(persistency_fit_fn=persistency_fn)

    # load image into ds9
    ds9.set_np2arr(rss.weighted_mean)

    if (display_filename is not None and os.path.isfile(display_filename)):
        # load the specified image into frame 2
        ds9.set("frame new")
        ds9.set("file %s" % (display_filename))

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
        kwargs=dict(rss=rss, ds9_queue=ds9_queue, save_plots=True),
    )
    plotter.daemon = True
    plotter.start()


    ds9_thread.join()
    plotter.join()

    # delete rss to clean up shared memory
    del rss

    print("all done!")
