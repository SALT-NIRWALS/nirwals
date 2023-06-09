#!/usr/bin/env python3

import matplotlib
print(matplotlib.get_backend())
# matplotlib.use('QT5Agg')
# matplotlib.use('GTK3Agg')
# matplotlib.use('WebAgg')
matplotlib.use('TkAgg')
# matplotlib.use('WebAgg')
print(matplotlib.get_backend())
import matplotlib.pyplot as plt
import numpy
import astropy.io.fits
import nirwals_reduce

x = numpy.linspace(-10,10,25)

plt.ion()
fig = plt.figure()
fig.show()

ax = fig.add_subplot(111)

while (True):

    slope = float(input())
    print(slope)

    y = x * slope
    fig.clf()
    ax = fig.add_subplot(111)
    ax.scatter(x,y)
    # fig.draw()
    fig.canvas.draw_idle()

