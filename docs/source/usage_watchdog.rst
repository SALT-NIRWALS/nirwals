.. _nirwals_watchdog:
**********************************************
On-the-fly reduction using *nirwals_watchdog*
**********************************************

The basic function of the ``nirwals_watchdog`` tool is to provide a on-the-fly
data reduction and inspection in near real-time for use at the telescope.
Towards this goal, this tool includes:

- reference pixel correction
- non-linearity correction
- conversion of raw reads to observed count rates in counts per second
- proper handling of saturated pixels
- automatic display of the resulting data file in ds9, with communication
  via the SAMP protocol.

Note that this tool does NOT include a stand-alone SAMP server (at least
not yet; i.e. for the automatic ds9 functionality to work you need an active
SAMP Hub, e.g. by running topcat).

.. udpate::

    The most recent version of ds9 has a built-in SAMP Hub, so no external
    standalone hub is required.

To run, one needs to specify a directory to monitor, a staging dir, and some
optional options that modify some of the internal algorithm::

    nirwals_watchdog incoming_data_dir --nonlinearity=/some/where/nonlinpoly.fits \
       --stage=/some/other/dir


Available options
==================

A full listing of all supported options can be obtained by running ``nirwals_reduce --help``::

    ...> nirwals_watchdog --help
    usage: nirwals_watchdog [-h] [--stage STAGING_DIR] [--nonlinearity NONLINEARITY_FN]
                            [--refpixel REF_PIXEL_MODE] [--test TEST] [--nowait] [--shmem] [--write] directory

    positional arguments:
      directory             name of directory to watch

    options:
      -h, --help            show this help message and exit
      --stage STAGING_DIR   staging directory
      --nonlinearity NONLINEARITY_FN
                            non-linearity correction coefficients (3-d FITS cube)
      --refpixel REF_PIXEL_MODE
                            reference pixels mode [default: NO]
      --test TEST           test mode; syntax: --test=delay:@filelist
      --nowait              do not wait for SAMP server
      --shmem               use shared memory to display files in ds9
      --write               write reduced frames (even if not needed)


Options mirror the naming convention and functionality of ``nirwals_reduce`` wherever possible.

- ``--stage`` staging directory

- ``--nonlinearity`` specify the non-linearity correction parameter file

- ``--refpixel`` reference pixels mode (use _blockyslope2_ for best performance)

- ``--test`` Run a test, simulating the arrival of new files in the monitored
  directory. Syntax is ``--test=delay:@filelist``, with delay giving a delay time
  between successive frames in seconds, and filelists specifying a file containing
  a list of "newly arrived" files, with one file per line (lines starting twith #
  are ignored)

  Run-times were tested on a modern laptop (i7 CPU, 32 GB RAM). Using a "full"
  reduction mode, including non-linearity takes approx 0.4 to 0.5 seconds per frame,
  from finding the newly arrived frame to end of writing the final result file.
  Given the minimum read time of the NIRWALS instrument of ~0.7 seconds this should allow
  monitoring incoming data in effectively real time (i.e. the previous frame is displayed
  before the next read is fully read out).