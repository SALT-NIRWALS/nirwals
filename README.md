# nirwals_reduce - instrumental detrending pipeline for SALT NIRWALS

# Components

The installation includes both the actual python package, as well as links to a small number of 
standalone tasks to run different functions of the code. These include

* `nirwals_reduce.py` -- run the science pipeline on one or more frames, taking raw data from individual reads to 
science-ready data, calibrated in counts per second, and including both uncertainties and several quality flags.

* `nirwals_fit_nonlinearity.py` -- generate non-linearity correction coefficients from a sequence of flat-field images

* ```nirwals_stalker.py``` -- On-the-fly reduction tool that monitors a specified directory, performs a limited but fast
data reduction, and displays the resulting results in a ds9 viewer in near-real time.

# Installation and dependencies

The nirwals package requires the following standalone python packages to run, most of which are likely already installed
on any common astronomy system:

* astropy
* matplotlib
* numpy
* scipy
* pandas 
* sklearn (scikit-learn)
* multiparlog -- parallel processing safe logging 
* pyvo -- VO capabilities, including the SAMP protocol used only for the on-the-fly reduction tool

All dependencies should be installed automatically as part of any installation (see below), but can also be installed
individually.

**Note:** depending on your local python environment setup, the installation may require admin privileges to install the
standalone tools. 

## installation via pip (recommended)

Nirwals pipeline is available via pip. To install, simply run

```
pip install nirwals
```

new pip releases are handled automatically via GitHub, so this provides access to the latest stable version.

## Installation directly from GitHub source

Alternatively, you can install the pipeline directly from the source on github. 

```
git clone https://github.com/SALT-NIRWALS/nirwals.git
cd nirwals
python3 setup.py install
```

Unlike the pip method above, this provides access to the bleeding edge of the code.


# How to run: `nirwals_reduce` data reduction pipeline

> nirwals_reduce [options] file.fits

### Available options

  `--maxfiles=N` specifies the maximum number of files to open for a given 
  up-the-ramp group. This is mostly to limit RAM usage. Default is no limit.

  `--nonlinearity=file.fits` 
  Apply non-linearity corrections to the reference-pixel/first-read subtracted 
  dataset. The reference file should be a file generated via the 
  rssnir_fit_nonlinearity tool to contain the pixel-level corrections in the 
  correct format

  `--flat=flat.fits`
  Specify a flatfield frame. Not implemented yet.

  `--dark=dark.fits`
  Subtract a dark-current correction from the entire input data cube. Use 
  _rssnir_makedark.py_ to generate the dark calibration frame.

  `--output=_suffix_` 
  When generating the output filename, the specified _suffix_ is inserted into the 
  input filename. Example: for input file _rss_test.fits_ the output filename would 
  be _rss_test.suffix.fits. Default is "reduced".

  `--refpixel` 
  Use the reference pixel in the first & last 4 rows and columns to 
  subtraced an instrumental pedestal level off all the input data. If not specified 
  the first read is considered to contain this zero-exposure offset. 

  `--dumps` Mostly used for debugging. When provided the tool also writes a number
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

### Caveats and limitations

- Not yet supported are fowler-reads of any kind, in particular when combined with 
  up the ramp sampling.
- Watch out when running on large numbers of up-the-ramp samples to avoid running out
  of memory (RAM). At this time the tool is optimized towards computing time at the 
  expense of memory demand. If in doubt or to begin use the _--maxfiles_ option to limit the number
  the number of open files and thus the memory footprint.


# How to run: On-the-fly reduction and monitoring

The basic function of the `nirwals_stalker.py` tool is to provide a on-the-fly data reduction and inspection in near 
real-time for use at the telescope. Towards this goal, this tool includes:

* reference pixel correction
* non-linearity correction 
* conversion of raw reads to observed count rates in counts per second
* proper handling of saturated pixels
* automatic display of the resulting data file in ds9, with communication via the SAMP protocol.

Note that this tool does NOT include a stand-alone SAMP server (at least not yet; i.e. for the automatic ds9 
functionality to work you need an active SAMP Hub, e.g. by running topcat).

To run, one needs to specify a directory to monitor, a staging dir, and some optional options that modify some of the 
internal algorithm.

```
nirwals_stalker.py incoming_data_dir --nonlinearity=/some/where/nonlinpoly.fits --stage=/some/other/dir
```

### Available options

Options mirror the naming convention and functionality of `nirwals_reduce` wherever possible.

* `--stage` staging directory
* `--nonlinearity` specify the non-linearity correction parameter file
* `--refpixel` reference pixels mode (use blockyslope2 for best performance)
* `--test` Run a test, simulating the arrival of new files in the monitored directory. Syntax is 
`--test=delay:@filelist`, with delay giving a delay time between successive frames in seconds, and filelists specifying 
a file containing a list of "newly arrived" files, with one file per line (lines starting twith # are ignored)

Run-times were tested on a modern laptop (i7 CPU, 32 GB RAM). Using a "full" reduction mode, including non-linearity 
takes approx 0.4 to 0.5 seconds per frame, from finding the newly arrived frame to end of writing the final result file. 
Given the minimum read time of the NIRWALS instrument of ~0.7 seconds this should allow monitoring incoming data in 
effectively real time (i.e. the previous frame is displayed before the next read is fully read out). 
