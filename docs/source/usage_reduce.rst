*************************************
Reducing data using *nirwals_reduce*
*************************************

The most basic call of nirwals involves specyfing just a single frame from
the read sequence to be reduced. One valid example would be::

    nirwals_reduce [options] file.fits

[options] are optional and specify certain operation parameters. For a full list of
available options you can invoke *nirwals_reduce* either without any parameters, or
by specifying the `--help` flag::

    nirwals_reduce --help

This will present the following condensed overview. A more detailed explanation
for each parameter is presented below::

    ...> nirwals_reduce --help
    usage: nirwals_reduce [-h] [--maxfiles MAX_NUMBER_FILES]
                          [--nonlinearity NONLINEARITY_FN] [--flat FLATFIELD_FN]
                          [--dark DARK_FN] [--output OUTPUT_POSTFIX]
                          [--persistency PERSISTENCY_MODE] [--saturation SATURATION]
                          [--dumps WRITE_DUMPS] [--debugpngs]
                          [--refpixel REF_PIXEL_MODE] [--flat4salt] [--report]
                          [--speedy] [--ncores N_CORES]
                          [--algorithm {rauscher2007,linreg,pairwise}]
                          files [files ...]

    positional arguments:
      files                 list of input filenames

    options:
      -h, --help            show this help message and exit
      --maxfiles MAX_NUMBER_FILES
                            limit number of files to load for processing
      --nonlinearity NONLINEARITY_FN
                            non-linearity correction coefficients (3-d FITS cube)
      --flat FLATFIELD_FN   calibration flatfield
      --dark DARK_FN        calibration dark
      --output OUTPUT_POSTFIX
                            addition to output filename
      --persistency PERSISTENCY_MODE
                            persistency mode
      --saturation SATURATION
                            saturation value/file
      --dumps WRITE_DUMPS   write intermediate process data [default: NO]
      --debugpngs           generate debug plots for all pixels with persistency [default: NO]
      --refpixel REF_PIXEL_MODE
                            reference pixels mode [default: NO]
      --flat4salt           write a flat, 1-extension FITS file for SALT
      --report              report ata provenance at end of processing
      --speedy              speed up processing by adaptively reducing sample reads
      --ncores N_CORES      number of CPU cores to use
      --algorithm {rauscher2007,linreg,pairwise}
                            number of CPU cores to use


Available options
=================

- ``maxfiles=N`` specifies the maximum number of files to open for a given up-the-ramp
  group. This is mostly to limit RAM usage. Default is no limit.

- ``--nonlinearity=file.fits`` Apply non-linearity corrections to the
  reference-pixel/first-read subtracted dataset. The reference file should be a file generated via the rssnir_fit_nonlinearity tool to contain the pixel-level corrections in the correct format

- ``--output=_suffix_`` When generating the output filename, the specified suffix
  is inserted into the input filename. Example: for input file rss_test.fits the
  output filename would be _rss_test.suffix.fits. Default is "reduced".

- ``--refpixel`` Use the reference pixel in the first & last 4 rows and columns
  to subtraced an instrumental pedestal level off all the input data. If not
  specified the first read is considered to contain this zero-exposure offset.

- ``--algorithm`` selects the algorithm to be used to estimate the observed signal rate
  from the individual up-the-ramp sequences.

- ``--dumps`` Mostly used for debugging. When provided the tool also writes a
  number of intermediate data products to disk that allow testing and verification.

- ``--report`` adds a condensed summary of all files and tools used during the reduction
  (that includes both software setup, input files, and specified calibration products).
  The resulting information is also stored as FITS table in the output file and can be
  read using the ``nirwals_provenance`` tool (see below for usage).

- ``--ncores=#`` allows to specify the number of parallel computing processes used during
  the processing. By default all available cores are used.

Additional options that are currently not fully implemented and/or operational, and/or
only to be used for debugging:

- ``--flat=flat.fits`` Specify a flatfield frame. Not implemented yet.

- ``--dark=dark.fits`` Subtract a dark-current correction from the entire input
  data cube. Use rssnir_makedark.py to generate the dark calibration frame. Currently
  not implemented since unstable dark-currents do not improve output data quality.


Example call::

    nirwals_reduce --refpixel=blockyslope2 --maxfiles=70 \
        SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.20.fits

output::

    rkotulla@legion:/work/rss/salt> ../rss_reduce/rss_reduce.py --refpixel \
        --maxfiles=70 SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.20.fits
    /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.20.fits
    /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.1.fits
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.2.fits
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.3.fits
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.4.fits
    ...
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.247.fits
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.248.fits
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.249.fits
     -- /work/rss/salt/SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.250.fits
    Limiting filelist to 70 frames
    (70, 2048, 2048)
    Applying non-linearity corrections
    No nonlinearity corrections loaded, skipping
    No linearized data found, using raw data instead
    No dark correction requested, skipping
    diff stack: (70, 2048, 2048)
    Identifying bad pixels
    Cleaning image cube
    calculating final image from stack
    Writing reduced results to 20220606_RN_URG_2reads_9dB.540.1.reduced.fits
    all done!


Caveats and limitations
==========================

- Not yet supported are fowler-reads of any kind, in particular when combined with up the ramp sampling.

- Watch out when running on large numbers of up-the-ramp samples to avoid running out of memory (RAM).
  At this time the tool is optimized towards computing time at the expense of memory demand. If in doubt or to
  begin use the ``--maxfiles`` option to limit the number the number of open files and thus the memory footprint.


Logging and more detailed progress updates
===========================================

