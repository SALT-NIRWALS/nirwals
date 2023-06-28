#!/usr/bin/env python3

import os
import sys
import numpy
import astropy.io.fits as pyfits

edge = 1

if __name__ == "__main__":

    sequence_file = sys.argv[1]
    output_file = sys.argv[2]

    filelist = []
    fileparts = sequence_file.split(".")
    hdulist = pyfits.open(sequence_file)
    n_reads = hdulist[0].header['NGROUPS']
    for read in range(1, n_reads+1):
        fn = "%s.%d.%s" % (".".join(fileparts[:-2]), read, fileparts[-1])
        print(fn)
        if (os.path.isfile(fn)):
            filelist.append(fn)
        else:
            break

    # print("\n".join([]+filelist))
    # print("output file: ", output_file)

    all_tops = []
    all_bottoms = []
    for fn in filelist:
        print(fn)
        hdulist = pyfits.open(fn)
        data = hdulist[0].data

        #
        # # first, combine left & right to subtract row-wise overscan level
        # _left = numpy.mean(data[:, edge:4], axis=1).reshape((-1, 1))
        # _right = numpy.mean(data[:, -4:-edge], axis=1).reshape((-1, 1))
        # row_wise = numpy.mean([_left, _right], axis=0)
        # if (debug):
        #     print(row_wise.shape)
        #
        # # plt.scatter(numpy.arange(row_wise.shape[0]), row_wise, s=1)
        # # plt.show()
        #
        # data_rowsub = data - row_wise
        # if (debug):
        #     pyfits.PrimaryHDU(data=data_rowsub).writeto("del__rowsub_%d.fits" % (dummycounter), overwrite=True)
        #
        # # now figure out the column situation
        top = data[edge:4, :]
        bottom = data[-4:-edge, :]

        mean_top = numpy.mean(top, axis=0)
        mean_bottom = numpy.mean(bottom, axis=0)

        # print(top.shape, mean_top.shape)

        all_tops.append(mean_top)
        all_bottoms.append(mean_bottom)
        # break

    all_tops = numpy.array(all_tops)
    all_bottoms = numpy.array(all_bottoms)

    print("Saving results to %s" % (output_file))
    pyfits.HDUList([
        pyfits.PrimaryHDU(),
        pyfits.ImageHDU(data=all_tops, name="TOPS"),
        pyfits.ImageHDU(data=all_bottoms, name="BOTTOMS"),
        pyfits.ImageHDU(data=(all_tops-all_bottoms), name="TOPS-BOTTOMS"),
    ]).writeto(output_file, overwrite=True)


        # ref_columns = numpy.vstack([top, bottom])
        # if (debug):
        #     print("ref-columns shape:", ref_columns.shape)
        # n_rows = 8 - 2 * edge
        # n_amps = 32
        # amp_size = ref_columns.shape[1] // n_amps
        #
