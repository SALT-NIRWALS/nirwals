#!/usr/bin/env python3

import sys
import os
import nirwals_reduce


if __name__ == "__main__":

    fn = sys.argv[1]

    try:
        saturation_fn = sys.argv[2]
    except:
        saturation_fn = 62000

    rss = rss_reduce.RSS(fn, saturation=saturation_fn)

    # rss.reduce(write_dumps=False)
    # rss.write_results()

    rss.load_all_files()
    # rss.subtract_first_read()

    rss.fit_nonlinearity(ref_frame_id=4, make_plot=False)

    # rss.plot_pixel_curve(818,1033)
    # rss.plot_pixel_curve(1700,555)
