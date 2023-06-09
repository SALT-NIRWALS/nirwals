# nirwals_reduce - instrumental detrending pipeline for SALT NIRWALS

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
    
### Example call:

```/work/rss/rss_reduce.py  --refpixel --maxfiles=70 SALT_data_RN_20220606/20220606_RN_URG_2reads_9dB.540.1.20.fits```

output:

```
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
```

## Caveats and limitations

- Not yet supported are fowler-reads of any kind, in particular when combined with 
  up the ramp sampling.
- Watch out when running on large numbers of up-the-ramp samples to avoid running out
  of memory (RAM). At this time the tool is optimized towards computing time at the 
  expense of memory demand. If in doubt or to begin use the _--maxfiles_ option to limit the number
  the number of open files and thus the memory footprint.