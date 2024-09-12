*******************************************
Creating non-linearity coefficient files
*******************************************

To generate nonlinearity files the pipeline comes with a special tool, ``nirwals_fit_nonlinearity``.
The tool is made to work virtually automated, with very few user inputs required. The most basic call  would be::

    nirwals_fit_nonlinearity [options] one_file_from_sequence.fits

[options] are optional and specify certain operation parameters. For a full list of
available options you can invoke ``nirwals_fit_nonlinearity`` either without any parameters, or
by specifying the `--help` flag::

    nirwals_fit_nonlinearity --help

This will present the following condensed overview. A more detailed explanation
for each parameter is presented below::

    usage: nirwals_fit_nonlinearity [-h] [--maxfiles MAX_NUMBER_FILES]
           [--nonlinearity NONLINEARITY_FN] [--saturation SATURATION] [--no_optimize]
           [--reflevel REFLEVEL] [--ncores N_CORES] [--refpixel REF_PIXEL_MODE]
           [--verify] [--pixels PIXELS] files

    positional arguments:
      files                 list of input filenames

    options:
      -h, --help            show this help message and exit
      --maxfiles MAX_NUMBER_FILES
                            limit number of files to load for processing
      --nonlinearity NONLINEARITY_FN
                            non-linearity correction coefficients (3-d FITS cube)
      --saturation SATURATION
                            saturation value/file
      --no_optimize         optimize saturation level as part of fitting
      --reflevel REFLEVEL   saturation value/file
      --ncores N_CORES      number of CPU cores to use for parallel fitting
      --refpixel REF_PIXEL_MODE
                            reference pixels mode [default: NO]
      --verify              verify results rather than fitting coefficients
      --pixels PIXELS       list of pixel coordinates or file with coordinates



Available options
=================

:kbd:`--maxfiles=N`
  specifies the maximum number of files to open for a given up-the-ramp
  group. This is mostly to limit RAM usage. Default is no limit.

:kbd:`--nonlinearity=file.fits`
  Determine the name of the __output__ nonlinearity correction file. This is made
  to be compatible with nirwals_reduce.

:kbd:`--saturation=some_level`
  Override the default saturation level used for the first initial fit. Note that
  nirwals_fit_nonlinearity automatically determines the saturation level for each
  pixel (unless the --no_optimize option is given, see below), so this can be
  safely left at the detault value

:kbd:`--no_optimize`
  Disable the algorithm that iteratively determines the saturation level of each pixel
  and along with it the full well depth.

:kbd:`--reflevel`
  Select a reference level used to find the illumination level (in counts/second)
  for each pixel. Default is 5000.

:kbd:`--ncores=#`
  allows to specify the number of parallel computing processes used during
  the processing. By default all available cores are used.

:kbd:`--refpixel`
  Use the reference pixel in the first & last 4 rows and columns
  to subtraced an instrumental pedestal level off all the input data. If not
  specified the first read is considered to contain this zero-exposure offset.
  Uses the same options as ``nirwals_reduce``

:kbd:`--verify`
  Used in conjunction with the ``--pixels`` option to plot the nonlinearity curve
  for select pixels.



Output
===========

The output of ``nirwals_fit_nonlinearity`` is a multi-extension FITS file containing
all information required to correct for detector nonlinearity. Extensions are as
follows:

* ``NONLINPOLY``

  | This 3-D image contains polynomial fitting coefficients for each pixel, ordered
  | from highest to lowest order term. The last two terms (linear term and constant
  | offset) are not used during the correction; the linear term is recalculated to be
  | identical to 1 so as not to affect noise and flatfield characteristics.

* ``FLAGS``

  Flag image containing a binary-coded representation of the outcome of the fitting
  procedure. If multiple flags apply to a given pixel then the FLAGS extension
  contains the OR-combined value. FLAGS value of 0 means everything worked as intended.

  Individual flags are:

  * NONLIN_FLAG_OK = 0x00

    no adverse findings, good value

  * NONLIN_FLAG_INSUFFICIENT_DATA = 0x01

    insufficient data for polynomial fit (less than 10 valid reads).
    This typically happens when pixels
    saturate too early, or if the --maxfiles option was set too low.

  * NONLIN_FLAG_NEGATIVE = 0x02

    When applying the nonlinearity fit to correct raw data, some corrected values
    are negative, indicating an illegal solution.

  * NONLIN_FLAG_BADSLOPE = 0x04

    When applying the nonlinearity fit to correct raw data, the relationship
    flux vs time encounters an area of negative slope, and thus can not be used
    to correct data.

  * NONLIN_FLAG_NEGATIVE_REFSLOPE = 0x08

    During the initial, unconstrained polynomial fit the linear term was negative.
    Likely cause is either bad data, or a saturation limit set too high (in which
    case all saturated data will be considered valid data, and the data might be
    dominated by a flat flux vs time relation, leading to a distorted fit).

  * NONLIN_FLAG_FITERROR = 0x10

    any other fit error

* ``SATURATION_LEVEL``

  raw intensity level for each pixel where the raw data is no longer linearizable.


* ``FULL_WELL_DEPTH_RAW``

  intensity above reference pixel level where data can no longer be linearized.


* ``FULL_WELL_DEPTH_CORRECTED``

  similar as FULL_WELL_DEPTH_RAW, but AFTER applying the nonlinearity correction.


* ``PRECISION_MEDIAN_CLIPPED``

  Precision is defined as non-linearity corrected flux divided by the perfect
  solution derived from the reflevel slope (without scatter or noise the precision
  would be identical to 1 at all intensities), after excluding outliers identified
  via iterative 3-sigma clipping.


* ``PRECISION_SIGMA_CLIPPED``

  Remaining scatter in the precision data, after sigma-clipping. Typical values
  should be <0.01, meaning that nonlinearity correction is able to linearize the
  data with better than 1% precision.


* ``PRECISION_SIGMA_FULL``

  Like PRECISION_SIGMA_CLIPPED, but before sigma clipping


* ``MAX_READ_RAW``

  Maximum raw intensity (before reference pixel correction) encountered during
  fitting. This value should be at or very close to saturation.


* ``MAX_READ_REFPIXCORR``

  Similar to MAX_READ_RAW, but AFTER reference pixel correction


* ``PROVENANCE``

  Data provenance, i.e. a full inventory of files read to generate this
  nonlinearity coefficient file.
