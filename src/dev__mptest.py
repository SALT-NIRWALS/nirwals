#!/usr/bin/env python3


import os
import sys

import astropy.io.fits as pyfits

import nirwals_reduce


def fitfct(rss, x, y):
    return rss.image_stack[1, y-1, x-1]

if __name__ == "__main__":

    fn = sys.argv[1]

    rss = rss_reduce.RSS(fn, max_number_files=10)
    rss.load_all_files()
    # rss.apply_nonlinearity_corrections()
    rss.image_stack -= rss.image_stack[0, :, :]

    results = rss.parallel_fitter(
        #xrange=[300, 400], yrange=None, #[200,275],
        execute_function=rss.fit_2component_persistency_plus_signal,
        is_in_class=True,
        return_dim=2,
    )

    print(results.shape)

    out_fn = sys.argv[2]
    print("Writing results to %s" % (out_fn))
    pyfits.PrimaryHDU(data=results).writeto(out_fn, overwrite=True)

