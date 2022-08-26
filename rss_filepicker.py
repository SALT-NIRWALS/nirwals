

import numpy
import os
import datetime
import astropy.io.fits as pyfits
import glob



class PreviousFilePicker( object ):

    def __init__(self, search_dir,autosave=True, update=False):
        self.search_dir = search_dir
        self.index_fn = os.path.join(search_dir, ".nirs_index.txt")
        self.index = {}

        self.load_index()
        self.inventory_directory(update=update)
        if (autosave):
            self.write_index()

    def load_index(self):
        if (not os.path.isfile(self.index_fn)):
            return
        with open(self.index_fn, "r") as f:
            lines = f.readlines()
            for line in lines:
                fn, timestamp_str = line.split()
                timestamp = datetime.datetime.fromisoformat(timestamp_str)
                self.index[fn] = timestamp
        pass

    def inventory_directory(self, update=False):
        filelist = glob.glob("%s/*.fits" % (self.search_dir))
        for fn in filelist:
            _, bn = os.path.split(fn)
            if (fn not in self.index or update):
                # this is a new file
                hdulist = pyfits.open(fn)
                timestamp = self.get_timestamp_from_hdr(hdulist[0].header)
                self.index[bn] = timestamp
        return len(self.index)


    def write_index(self):
        with open(self.index_fn, "w") as f:
            for fn in self.index:
                print("%s %s" % (fn, self.index[fn].isoformat()), file=f)

    def get_timestamp_from_hdr(self, hdr):
        date_obs = hdr['DATE-OBS']
        # print(date_obs, type(date_obs))
        items = date_obs.split()
        iso_fmt = "%sT%s+0%s" % (items[0], items[1], items[3])
        dt = datetime.datetime.fromisoformat(iso_fmt)
        return dt

    def find_closest(self, hdr):
        timestamp = self.get_timestamp_from_hdr(hdr)

        closest_dt = None
        closest_fn = None
        for fn in self.index:
            datetime = self.index[fn]
            delta_t = datetime - timestamp
            delta_seconds = delta_t.total_seconds()
            # print(delta_t, fn)

            if (closest_dt is None and delta_seconds < 0):
                # this is the first encountered prior frame
                closest_dt = numpy.fabs(delta_t.total_seconds())
                closest_fn = fn
            elif (closest_dt is None):
                # not prior
                continue
            elif (delta_seconds < 0 and numpy.fabs(delta_seconds) < closest_dt):
                # new closest match
                print("updating closest match: % 10.3f :: %s" % (delta_seconds, fn))
                closest_dt = numpy.fabs(delta_seconds)
                closest_fn = fn

        if (closest_fn is not None):
            closest_fn = os.path.abspath(os.path.join(self.search_dir, closest_fn))
        return closest_fn, closest_dt

