#!/usr/bin/env python3

import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import itertools
import os
from matplotlib.backends.backend_pdf import PdfPages

import sys
import os

def make_textpage(txt, pdf, fontsize=24):

    firstPage = plt.figure()# %figsize=(11.69,8.27))
    firstPage.clf()
    # txt = 'This is the title page'
    firstPage.text(0.5,0.5,txt, transform=firstPage.transFigure, size=fontsize, ha="center")
    pdf.savefig(firstPage)
    plt.close()




if __name__ == '__main__':

    input_fn = sys.argv[1]
    items = input_fn.split(".")
    items[-2] = '%d'
    bn = ".".join(items)
    # bn = "../SALT_data_2024/0607/raw/N202406070006.6.2.%d.fits"

    pdf_fn = sys.argv[2]
    if (pdf_fn == 'auto'):
        _, x = os.path.split(bn)
        items = x.split(".")
        pdf_fn = "%s_%s.pdf" % (items[0], items[2])
    print("Preparing output PDF: %s" % (pdf_fn))
    if (os.path.isfile(pdf_fn)):
        print("Output file already exists")
        sys.exit(-1)

    output_hdu_fn = "%s_%s_refpixelstats.fits" % (items[0], items[2])

    pdf = PdfPages(pdf_fn)

    # title page with input filename
    _, base = os.path.split(bn)
    hdulist = pyfits.open(input_fn)
    try:
        objectname = hdulist[0].header['OBJECT']
    except:
        objectname = "unknown"
    make_textpage("%s\nOBJECT: %s" % (base, objectname), pdf)


    #
    # Read all data
    #
    print("Reading data")
    refpixel_stack = []
    for read in range(1, 824):
        fn = bn % read
        # print(fn, os.path.isfile(fn))

        if (not os.path.isfile(fn)):
            continue

        print(">", end='')
        if ((read % 100) == 0):
            print("\n", end='')
        elif ((read % 10) == 0):
            print(" ", end='')

        hdulist = pyfits.open(fn)
        data = hdulist[0].data.astype(float)
        _top = data[:4, :]
        _bottom = data[-4:, :]
        refpix = numpy.vstack([_top, _bottom])
        # print(refpix.shape)

        refpixel_stack.append(refpix)
    print(" done!")

    refpixel_stack = numpy.array(refpixel_stack)

    # work out number of reads
    n_reads = refpixel_stack.shape[0]
    reads = numpy.arange(n_reads)
    print(n_reads)

    # take out average level of each pixel
    per_pixel_average = numpy.mean(refpixel_stack, axis=0)
    #print(per_pixel_average.shape)
    pixel_avg_subtracted = refpixel_stack-per_pixel_average
    #print(pixel_avg_subtracted.shape)


    #
    # Stats for each amplifier to identify average trends
    #
    n_amps = 32

    amp_resorted = numpy.swapaxes(pixel_avg_subtracted, 1, 2)
    print("amp resorted:", amp_resorted.shape)

    per_amp = amp_resorted.reshape((n_reads, n_amps, -1))
    print("per amp:", per_amp.shape)
    per_amp_average = numpy.median(per_amp, axis=2)
    print("per amp average: ", per_amp_average.shape)

    per_amp_fullres = numpy.repeat(per_amp_average, repeats=64, axis=1).reshape((n_reads, 1, -1))
    print("per amp fullres:", per_amp_fullres.shape)


    amp_avg_subtracted = pixel_avg_subtracted - per_amp_fullres
    print("amp avg subtracted:", amp_avg_subtracted.shape)


    # Do some smoothing over time/reads
    window = 20
    padded = numpy.pad(per_amp_average, pad_width=((window,window), (0,0)), mode='constant', constant_values=numpy.nan)
    print(per_amp_average.shape, padded.shape)

    # print(n_reads)
    padded_3d = numpy.empty((2*window+1, per_amp_average.shape[0], per_amp_average.shape[1]))
    print(padded_3d.shape)
    for i, offset in enumerate(range(-window, window+1)):
        # print(i, offset)
        padded_3d[i,:,:] = padded[window+offset:n_reads+window+offset, :]

    stats = numpy.nanpercentile(padded_3d, [16,50,84], axis=0)
    print(stats.shape)

    peramp_median = stats[1,:,:]
    peramp_sigma  = 0.5*(stats[2,:,:] - stats[0,:,:])



    peramp_global_mean = numpy.nanmean(per_amp_average, axis=0)
    peramp_global_median = numpy.nanmedian(per_amp_average, axis=0)
    peramp_global_std = numpy.nanstd(per_amp_average, axis=0)
    print(peramp_global_mean)
    print(peramp_global_median)
    print(peramp_global_std)


    #
    # Plot raw and average data for each amplifier
    #
    for amp in range(32):
        print(amp+1, end=' ')
        fig, ax = plt.subplots()
        ax.scatter(reads, per_amp_average[:,amp], s=0.3)
        ax.set_title("%s -- AMP %d" % (bn, amp+1))
        ax.set_xlabel("Read-#")
        ax.set_ylabel("amp-average mean level")

        ax.plot(reads, peramp_median[:,amp], lw=3, alpha=0.6)
        ax.plot(reads, peramp_median[:,amp]-peramp_sigma[:,amp], lw=1, alpha=0.6)
        ax.plot(reads, peramp_median[:,amp]+peramp_sigma[:,amp], lw=1, alpha=0.6)

        ax.plot(reads, peramp_sigma[:,amp], lw=1, alpha=1, c='red')

        pdf.savefig(fig)
        plt.close(fig)
        # break
    print("   :: all done!")



    #
    # Show a couple of example pixels
    #
    make_textpage("Example pixels", pdf)
    for amp in range(32)[::4]:
        fig, ax = plt.subplots()
        ax.scatter(reads, per_amp_average[:, amp], s=0.5)
        ax.scatter(reads, pixel_avg_subtracted[:,5,amp*64+27], s=0.5)

        ax.set_title("AMP %d" % (amp+1))

        aas = amp_avg_subtracted[:,5,amp*64+28]
        _stats = numpy.nanpercentile(aas, [16,50,84])
        _median = _stats[1]
        _sigma = 0.5 * (_stats[2]-_stats[0])
        print(_median, _sigma)
        ax.scatter(reads, aas, s=0.5)
        ax.axhline(y=0)

        pdf.savefig(fig)
        plt.close(fig)



    # identify and show outliers with more noise than there should be
    _stats = numpy.nanpercentile(amp_avg_subtracted, [16,50,84], axis=0)
    print(_stats.shape)

    _median = _stats[1]
    _sigma = 0.5*(_stats[2]-_stats[0])
    print(_median.shape)

    # find outliers
    outlier_limit = 25
    iy,ix = numpy.indices(_median.shape)
    iy[4:,:] += 2040

    large_sigma = _sigma > outlier_limit
    bad_x = ix[large_sigma]
    bad_y = iy[large_sigma]
    outlier_xy = numpy.vstack([bad_x, bad_y]).T
    # print(outlier_xy)
    print("Found %d outlier" % (outlier_xy.shape[0]))

    make_textpage("noise in each pixel\nafter removing amplifier average", pdf)

    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, ax = plt.subplots(figsize=(20, 5))
    x = numpy.arange(refpixel_stack.shape[2])
    for y in range(8):
        typ = numpy.nanmedian(_sigma[y, :])
        print(y, typ)
        # ax.scatter(x, _sigma[y,:]+y*10, s=0.2)
        label = "Y=%d" % (y + 1 if y < 4 else y + 2041)

        is_outlier = (_sigma[y, :] > outlier_limit)
        ax.scatter(x[~is_outlier], _sigma[y, :][~is_outlier], c=prop_cycle[y], s=0.2, label=label)
        ax.scatter(x[is_outlier], _sigma[y, :][is_outlier], c=prop_cycle[y], s=25)
        # ax.axhline(y=typ+y*10)
        # print(_median[y,:20])
    ax.axhline(y=outlier_limit, ls='--')
    ax.set_xlabel("Read")
    ax.set_ylabel("noise")
    ax.set_xlim(-20, 2250)

    ax.legend(markerscale=15)
    pdf.savefig(fig)

    fig, ax = plt.subplots(figsize=(20, 5))
    x = numpy.arange(refpixel_stack.shape[2])
    for y in range(8):
        typ = numpy.nanmedian(_sigma[y, :])
        print(y, typ)
        label = "Y=%d" % (y + 1 if y < 4 else y + 2041)
        ax.scatter(x, _sigma[y, :] + (7 - y) * 20, s=0.3, label=label)

    ax.set_xlabel("Read")
    ax.set_ylabel("noise + offset [counts]")
    ax.set_xlim(-20, 2250)
    ax.set_ylim((0, 180))
    ax.legend(markerscale=15)
    pdf.savefig(fig)



    make_textpage("Selected outliers based on noise", pdf)

    for x,y in outlier_xy[:20]:
        # print(x,y)
        _y = y if y <= 4 else y-2040
        amp = int(numpy.ceil(x/64.))
        fig,ax = plt.subplots(figsize=(10,4))
        ax.scatter(reads, pixel_avg_subtracted[:,_y,x], s=0.5, c='blue', label='pixel avg sub')
        # ax.scatter(reads, per_amp_fullres[:, 0, x], s=0.5, label='per amp average')
        ax.plot(reads, peramp_median[:,amp], c='blue', alpha=0.4, label="amp average")
        ax.scatter(reads, amp_avg_subtracted[:,_y,x], s=0.5, c='orange', label='amp avg sub')
        ax.set_title("x=%d  y=%d" % (x,y))
        ax.legend(markerscale=10)

        pdf.savefig(fig)


    # close pdf
    pdf.close()

    print("Saving all results")
    raw_hdu = pyfits.ImageHDU(data=refpixel_stack, name='REFPIXEL_STACK')
    per_pixel_average_hdu = pyfits.ImageHDU(data=per_pixel_average, name='PER_PIXEL_AVERAGE')
    pixel_avg_subtracted_hdu = pyfits.ImageHDU(data=pixel_avg_subtracted, name='PER_PIXEL_AVG_SUBTRACTED')
    per_amp_average_hdu = pyfits.ImageHDU(data=per_amp_average, name='PER_AMP_AVERAGE')
    per_amp_fullres_hdu = pyfits.ImageHDU(data=per_amp_fullres, name='PER_AMP_FULLRES')

    stats_hdu = pyfits.ImageHDU(data=stats, name='STATS')
    peramp_median_hdu = pyfits.ImageHDU(data=peramp_median, name='PER_AMP_MEDIAN')
    peramp_sigma_hdu = pyfits.ImageHDU(data=peramp_sigma, name='PER_AMP_SIGMA')
    outliers_hdu = pyfits.ImageHDU(data=outlier_xy, name='OUTLIERS_XY')
    outliers_hdu.header['LIMIT'] = outlier_limit

    prim_hdu = pyfits.PrimaryHDU(header=hdulist[0].header)
    prim_hdu.header['INPUT_FN'] = input_fn
    results_hdu = pyfits.HDUList([
        prim_hdu,
        raw_hdu,
        per_pixel_average_hdu, pixel_avg_subtracted_hdu,
        per_amp_average_hdu, per_amp_fullres_hdu,
        stats_hdu, peramp_median_hdu, peramp_sigma_hdu,
        outliers_hdu
    ])
    results_hdu.writeto(output_hdu_fn, overwrite=True)

    print("all done")

