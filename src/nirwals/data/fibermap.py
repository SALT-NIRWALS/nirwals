#!/usr/bin/env python3
import astropy.table
import pandas
import astropy.io.fits as pyfits
import sys
import os
import pandas




def get_path():
    _file = os.path.abspath(__file__)
    _dir, _ = os.path.split(_file)
    return _dir

def get_fibermap():
    csv_fn = os.path.join( get_path(), "fiber_map_20221018.csv")
    return csv2fits(csv_fn)

def csv2fits(csv_fn):
    csv = pandas.read_csv(csv_fn)
    table = astropy.table.Table.from_pandas(csv)
    tbhdu = pyfits.BinTableHDU(table, name="FIBERMAP")
    return tbhdu






if __name__ == "__main__":

    opt = sys.argv[1]

    if (opt == "csv2fits"):
        csv_fn = sys.argv[2]
        fits_fn = sys.argv[3]
        csv2fits(csv_fn, fits_fn)

    else:
        print("Not sure what to do here")
