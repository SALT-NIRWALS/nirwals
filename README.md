# rss_reduce - instrumental detrending pipeline for SALT RSS-NIR

## How to run

> rss_reduce [options] file.fits

### Available options

  **--maxfiles=N** specifies the maximum number of files to open for a given 
  up-the-ramp group. This is mostly to limit RAM usage. Default is no limit.

  **--nonlinearity=file.fits** 
  Apply non-linearity corrections to the reference-pixel/first-read subtracted 
  dataset. The reference file should be a file generated via the 
  rssnir_fit_nonlinearity tool to contain the pixel-level corrections in the 
  correct format

  **--flat=flat.fits**
  Specify a flatfield frame. Not implemented yet.

  **--dark=dark.fits**
  Subtract a dark-current correction from the entire input data cube. Use 
  _rssnir_makedark.py_ to generate the dark calibration frame.

  **--output=_suffix_** 
  When generating the output filename, the specified _suffix_ is inserted into the 
  input filename. Example: for input file _rss_test.fits_ the output filename would 
  be _rss_test.suffix.fits. Default is "reduced".

  **--refpixel** 
  Use the reference pixel in the first & last 4 rows and columns to 
  subtraced an instrumental pedestal level off all the input data. If not specified 
  the first read is considered to contain this zero-exposure offset. 

  **--dumps** Mostly used for debugging. When provided the tool also writes a number
  of intermediate data products to disk that allow testing and verification.
    