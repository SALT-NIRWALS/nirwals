#!/usr/bin/env python3


import numpy
import astropy.io.fits as pyfits
import os
import sys




if __name__ == "__main__":

    out_fn = sys.argv[1]


    in_filelist = sys.argv[2:]

    raw_cubes = []
    for fn in in_filelist:
        print("Reading input from %s" % (fn))
        hdulist = pyfits.open(fn)
        _cube = hdulist[0].data
        raw_cubes.append(_cube)

    n_files = len(raw_cubes)

    n_reads = [cube.shape[0] for cube in raw_cubes]
    print("#reads: %s" % (str(n_reads)))
    max_reads = numpy.max(numpy.array(n_reads))
    print("working with a max of %d reads" % (max_reads))

    linechunks = 1
    n_chunks = int(numpy.ceil(raw_cubes[0].shape[1] / float(linechunks)))
    print("We'll split the work into %d chunks of %d lines each" % (n_chunks, linechunks))

    nx, ny = 2048, 2048
    output_median = numpy.full((max_reads, ny, nx), fill_value=numpy.NaN)
    output_sigma = numpy.full((max_reads, ny, nx), fill_value=numpy.NaN)

    # max_reads = 15
    for r in range(max_reads):

        readbuffer = []
        for frame in range(n_files):
            if (r < n_reads[frame]):
                readbuffer.append(raw_cubes[frame][r,:,:])
        print("frame %d: n-reads: %d" % (r, len(readbuffer)))

        readbuffer = numpy.array(readbuffer)
        print("readbuffer shape:", readbuffer.shape)

        good_data = numpy.isfinite(readbuffer)

        # no iterations for now
        # for iteration in range(1):
        #     print(r, good_data.shape)
        #     _stats = numpy.nanpercentile(readbuffer, [16,50,84], axis=0)
        #     print(len(_stats), _stats[1].shape)
        #     # break

        # median_level = readbuffer[0, :, :]
        median_level = numpy.nanmean(readbuffer, axis=0)
        sigma = numpy.sqrt(numpy.nanvar(readbuffer, axis=0))

        #_stats = numpy.nanpercentile(readbuffer, [16,50,84], axis=0)
        #median_level = _stats[1]
        #sigma = 0.5 * (_stats[2] - _stats[0])

        output_median[r] = median_level
        output_sigma[r] = sigma

        print()

    pyfits.PrimaryHDU(data=output_median).writeto("cubecombine__mean.fits", overwrite=True)
    pyfits.PrimaryHDU(data=output_sigma).writeto("cubecombine__sigma.fits", overwrite=True)

