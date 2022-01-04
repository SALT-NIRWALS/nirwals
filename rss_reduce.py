#!/usr/bin/env python3
import multiprocessing
import os
import sys
import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import itertools
import multiprocessing

from astropy import log
log.setLevel('ERROR')

import astropy
print(astropy.__path__)

# def fit_nonlinearity_sequence(pinit, args):
#     """
#     Mostly a wrapper around scipy.optimize.leastsq.
#
#     """
#
#     def fit_fct(p, x):
# #         y = numpy.zeros(x.shape)
# #         for i in range(p.shape[0]):
# #             y += p[i] * x**(i+1)
# #         return y
#     def err_fct(p,x,y,err, fitrange_x, fitrange_y):
#         yfit = fit_fct(p,x)
#         in_fit_range = numpy.isfinite(x) & numpy.isfinite(y)
#         if (fitrange_x is not None):
#             in_fit_range = in_fit_range & (x >= fitrange_x[0]) & (x < fitrange_x[1])
#         if (fitrange_y is not None):
#             in_fit_range = in_fit_range & (y >= fitrange_y[0]) & (y < fitrange_y[1])
#         if (err is None):
#             return ((y-yfit))[in_fit_range]
#         return ((y-yfit)/err)[in_fit_range]
#
#     # (medlevel, exptime, None, intensity_range, exptime_range) = args
#
#     fit = scipy.optimize.leastsq(err_fct, pinit, args=args, full_output=1)
#
#     pfit = fit[0]
#
#     uncert = numpy.sqrt(numpy.diag(fit[1]))
#     print(pfit, uncert)
#     return pfit, uncert
#


class RSS(object):

    def __init__(self, fn):
        self.fn = fn
        self.filelist = []

        self.image_stack_initialized = False
        self.first_read_subtracted = False
        self.first_read = None

        self.nonlin_fn = None
        self.nonlinearity_cube = None

        self.get_full_filelist()

    def get_full_filelist(self):
        # get basedir
        fullpath = os.path.abspath(self.fn)
        print(fullpath)
        self.basedir, _fn = os.path.split(fullpath)
        self.filebase = ".".join(_fn.split(".")[:-2])
        # print(self.basedir, filebase)

        for _read in range(1,100):
            filename = "%s.%d.fits" % (self.filebase, _read)
            full_filename = os.path.join(self.basedir, filename)
            if (os.path.isfile(full_filename)):
                self.filelist.append(full_filename)
            else:
                break

            # print(full_filename, os.path.isfile(full_filename))

        print("\n -- ".join(self.filelist))

        return

    def add_file(self, filename):
        return

    def load_all_files(self):

        if (self.image_stack_initialized):
            print("stack already initialized, skipping repeat try")
            return

        self._image_stack = []

        # open all frames
        for fn in self.filelist[:50]:
            hdulist = pyfits.open(fn)
            # hdulist.info()
            imgdata = hdulist[0].data
            self._image_stack.append(imgdata)
            # break

        # calculate the initial image stack
        self.image_stack = numpy.array(self._image_stack, dtype=numpy.float32)
        print(self.image_stack.shape)

        self.image_stack_initialized = True

    def reduce(self, write_dumps=False):

        self.load_all_files()

        # apply any necessary corrections for nonlinearity and other things
        self.subtract_first_read()

        # calculate differential stack
        linearized = self.apply_nonlinearity_corrections()
        if (linearized is None):
            print("No linearized data found, using raw data instead")
            linearized = self.image_stack

        #self.image_stack
        self.differential_stack = numpy.diff(linearized, axis=0)
        print(self.differential_stack.shape)

        # mask out all saturated and/or otherwise bad samples
        bad_data = (self.image_stack > 60000)
        self.clean_stack = self.differential_stack.copy()
        self.clean_stack[bad_data[1:]] = numpy.NaN

        # calculate a average countrate image
        image = numpy.nanmean(self.clean_stack, axis=0)
        print(image.shape)

        ratios = linearized / linearized[3:4, :, :]

        if (write_dumps):
            pyfits.PrimaryHDU(data=self.image_stack+self.first_read).writeto("stack_raw.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.image_stack).writeto("stack_zerosub.fits", overwrite=True)
            pyfits.PrimaryHDU(data=linearized).writeto("stack_linearized.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.differential_stack).writeto("stack_diff.fits", overwrite=True)
            pyfits.PrimaryHDU(data=self.clean_stack).writeto("stack_clean.fits", overwrite=True)
            pyfits.PrimaryHDU(data=ratios).writeto("stack_ratios.fits", overwrite=True)
            pyfits.PrimaryHDU(data=image).writeto("final_image.fits", overwrite=True)

        return

    def subtract_first_read(self):
        if (self.first_read_subtracted):
            print("First read already subtracted, skipping")
            return

        if (self.first_read is None):
            self.first_read = self.image_stack[0].copy()

        self.image_stack -= self.first_read
        self.first_read_subtracted = True


    def _nonlinearity_fit_fct(self, p, x):
        y = numpy.zeros(x.shape)
        for i in range(p.shape[0]):
            y += p[i] * x**(i+1)
        return y
    def _nonlinearity_fit_err_fct(self, p, x, y):
        yfit = self._nonlinearity_fit_fct(p, x)
        err = numpy.sqrt(y + 10 ** 2)
        return ((y - yfit) / err)

    def _fit_nonlinearity_pixel(self, _x, _y):
        # print("fitting", _x.shape, _y.shape)
        # return 1

        # define two sub-routines we need for the fitting


        # prepare the data and initial fitting parameters
        # print(ix, iy)
        # _x = masked[:, iy, ix]
        # _y = linearized_intensity[:, iy, ix]
        good4fit = numpy.isfinite(_x) & numpy.isfinite(_y)
        # print(_x)
        # print(_y)
        x = _x[good4fit]
        y = _y[good4fit]
        if (numpy.sum(good4fit) < 5):
            # if there's no good data we can't do any fitting
            return numpy.array([1., 0., 0.]) # assume perfect linearity

        pinit = [1., 0., 0.]
        fit = scipy.optimize.leastsq(
            func=self._nonlinearity_fit_err_fct, x0=pinit,
            args=(x, y),
            full_output=1
        )
        pfit = fit[0]

        return pfit

    def fit_nonlinearity(self, ref_frame_id=10, max_linear=50000, make_plot=False):

        self.subtract_first_read()
        if (self.first_read_subtracted):
            bad_data = (self.image_stack + self.first_read) > max_linear
        else:
            bad_data = self.image_stack > max_linear

        # normalize each image with the N-th frame to take out a linear slope
        # normalized = self.image_stack / (self.image_stack[ref_frame_id] / ref_frame_id)

        # mask out all pixels above a saturation level
        # normalized[bad_data] = numpy.NaN

        masked = self.image_stack.copy()
        masked[bad_data] = numpy.NaN

        linearized_intensity = numpy.arange(masked.shape[0]).reshape((-1,1,1)) * (self.image_stack[ref_frame_id:ref_frame_id+1] / ref_frame_id)
        print(linearized_intensity.shape)
        pyfits.PrimaryHDU(data=linearized_intensity).writeto("linearized.fits", overwrite=True)


        # now fit a each pixel

        parallel = False

        if (parallel):
            _iy,_ix = numpy.indices((masked.shape[1], masked.shape[2]))
            _ixy = numpy.dstack([_ix,_iy]).reshape((-1,2))
            print("ixy shape:", _ixy.shape)
            _ixy = _ixy[:25]
            # print(_ixy[:25])


            pool = multiprocessing.Pool(processes=2) #multiprocessing.cpu_count())

            __ixy = list(zip(_ixy[:,0], _ixy[:,1]))
            _masked = [masked[iy,ix] for (ix,iy) in __ixy]
            _linint = [linearized_intensity[iy,ix] for (ix,iy) in __ixy]
            # print(masked)
            # print(list(zip(_ixy[:,0], _ixy[:,1])))

            print("masked")
            print(_masked)
            print(len(_masked))
            print(len(list(_masked)))

            print("linint")
            print(_linint)
            print(len(list(_linint)))

            print("it")
            it = zip(_masked, _linint)
            print(it)
            print(len(list(it)))

            results_parallel = pool.starmap(
                self._fit_nonlinearity_pixel,
                iterable=zip(_masked, _linint),
            )
            pool.close()
            pool.join()

            # [masked[:, iy, ix], linearized_intensity[:, iy, ix] for [ix,iy] in _ixy],
            #     # iterable=zip(itertools.repeat(masked),
            #     #              itertools.repeat(linearized_intensity),
            #     #              _ixy[:, 0], _ixy[:, 1]),

            print("results_parallel=\n", results_parallel)
        else:
            nonlinearity_fits_3d = numpy.zeros((3, masked.shape[1], masked.shape[2]))
            nonlinearity_fits_inverse = numpy.zeros((3, masked.shape[1], masked.shape[2]))
            for (ix,iy) in itertools.product(range(masked.shape[1]),
                                             range(masked.shape[2])):
                if (iy == 0):
                    sys.stdout.write("\rWorking on column % 4d" % (ix))
                    sys.stdout.flush()

                    # print(ix, iy)
                # if ((ix % 100) == 0):
                #

                # sys.stdout.write("ix=% 4d  iy=% 4d\r" % (ix,iy))
                # sys.stdout.flush()
                # print(ix, iy)

                _x = masked[:, iy, ix]
                _y = linearized_intensity[:, iy, ix]
                pfit = self._fit_nonlinearity_pixel(_x, _y)
                if (pfit is not None):
                    nonlinearity_fits_3d[:, iy:iy+4, ix:ix+4] = pfit.reshape((-1,1,1))

                pfit_inverse = self._fit_nonlinearity_pixel(_y, _x)
                if (pfit is not None):
                    nonlinearity_fits_inverse[:, iy:iy+4, ix:ix+4] = pfit_inverse.reshape((-1,1,1))

                if (make_plot):
                    fig = plt.figure()
                    ax = fig.add_subplot()
                    ax.scatter(_x, _y)
                    ax.plot(_x, self._nonlinearity_fit_fct(pfit, _x))
                    fig.suptitle("y = %+.4g*x %+.4g*x^2 %+.4g*x^3" % (pfit[0], pfit[1], pfit[2]))
                    fig.savefig("nonlin__x%04d__y%04d.png" % (ix, iy))
                    plt.close(fig)

            #break

        # return

        pyfits.PrimaryHDU(data=nonlinearity_fits_3d).writeto("nonlin3d.fits", overwrite=True)
        pyfits.PrimaryHDU(data=nonlinearity_fits_inverse).writeto("nonlin_inverse.fits", overwrite=True)
        return

    def read_nonlinearity_corrections(self, nonlin_fn):

        if (os.path.isfile(nonlin_fn)):
            try:
                hdulist = pyfits.open(nonlin_fn)
                print("Reading nonlinearity corrections from %s" % (nonlin_fn))
                nonlinearity_cube = hdulist[0].data
                print("CORR shape:", nonlinearity_cube.shape)
            except:
                return False

        self.nonlin_fn = nonlin_fn
        self.nonlinearity_cube = nonlinearity_cube

    def apply_nonlinearity_corrections(self):

        if (self.nonlinearity_cube is None):
            print("No nonlinearity corrections loaded, skipping")
            return None

        self.subtract_first_read()

        img_cube = self.image_stack
        print("NONLIN: data=%s   corr=%s" % (str(img_cube.shape), str(self.nonlinearity_cube.shape)))

        linearized_cube = \
            self.nonlinearity_cube[0:1, :, :] * img_cube + \
            self.nonlinearity_cube[1:2, :, :] * numpy.power(img_cube, 2) + \
            self.nonlinearity_cube[2:3, :, :] * numpy.power(img_cube, 3)

        return linearized_cube


    def write_results(self, fn=None):
        return

    def plot_pixel_curve(self, x, y):

        self.subtract_first_read()
        counts = self.image_stack[:, y-1, x-1]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(numpy.arange(counts.shape[0]), counts)
        fig.savefig("pixelcurve_x%04d_y%04d.png" % (x,y))



if __name__ == "__main__":

    fn = sys.argv[1]
    rss = RSS(fn)
    rss.read_nonlinearity_corrections("nonlin_inverse.fits")
    # rss.read_nonlinearity_corrections("nonlin3d.fits")
    rss.reduce(write_dumps=True)
    rss.write_results()

    rss.plot_pixel_curve(818,1033)
    rss.plot_pixel_curve(1700,555)


    print("all done!")