*********************
Data Products
*********************

.. dataproducts:

Generated output files
**************************

Depending on user configuration, `nirwals_reduce` generates multiple output files. These are typically named
after the filename of the input file, followed by a standardized prefix to identify the data product:

* ***_reduced.fits** These are the reduced files, and are the only files generated in any case. See below
  for details on FITS extensions and data units

* ***__stack_raw.fits** are full-cube output files containing all input reads combined into a single data
  cube. The data in this file is straight from the detector, with no calibrations applied. Pixels brighter
  than the configured saturation level are masked as NaNs.

* ***__stack_refpixcorr.fits** Similiar to the _stack_raw.fits, but after applying the per-read reference
  pixel correction.

* ***__stack_linearized.fits** As above, but with reference pixel AND non-linearity correction applied.
  This is the final stage of data preparation used for the up-the-ramp fitting.


Output file extensions
*************************

The final reduced files contain a number of output extensions that contain results from different methods or
metadata describing and/or derived from and as part of the data processing.

* **PRIMARY**
  The primary FITS header does not contain any image data, but rather all headers relevant to describe the
  input and output data. Most of the FITS headers are taken from the selected reference exposure (typically
  the first in the read sequence), but also contains a number of keywords summarizing the entire sequence.
  See below for more details.

* **SCI**:
  This is the final calibrated slope derived from the up-the-ramp fitting.

  **UNITS** Data is given in data numbers (counts) per second, with NO gain correction.

* **MEDIAN**:
  tbd

* **NOISE**
  tbd

* **NPAIRS**

* **MAX_T_EXP**
  Per-pixel maximum integration before exceeding the saturation level.

* **PROVENANCE**
  Unlike all other extensions, this extension is a FITS table, with key/value pairs describing many aspects
  relevant for reproducing the output file. This includes:

  * List of input files

  * Calibration products used (e.g. non-linearity corrections)

  * System configuration

  * Package version information for the nirwals package itself and other important packages
    (such as numpy & astropy)





Important and/or derived FITS headers
*****************************************
