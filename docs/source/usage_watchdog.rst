.. _nirwals_watchdog:

**********************************************
On-the-fly reduction using *nirwals_watchdog*
**********************************************

The basic function of the :command:`nirwals_watchdog` tool is to provide a on-the-fly
data reduction and inspection in near real-time for use at the telescope.
Towards this goal, this tool includes:

- reference pixel correction
- non-linearity correction
- conversion of raw reads to observed count rates in counts per second
- proper handling of saturated pixels
- automatic display of the resulting data file in ds9, with communication
  via the SAMP protocol.


.. Note::

    Note that this tool does NOT include a stand-alone SAMP server (at least
    not yet; i.e. for the automatic ds9 functionality to work you need an active
    SAMP Hub, e.g. by running topcat).

    **UPDATE:** The most recent version of ds9 has a built-in SAMP Hub, so no external
    standalone hub is required.

.. Tip::

    By default, each frame is reduced, written to disk, and then a message is sent to ds9 to update the frame and
    display the new image. This, therefore, likely requires some disk I/O, although in practice the operating system
    should cache the file, and this step not actually be slowed down by actually writing to and reading from a physical
    disk.

    Alternatively, the tool can interact with ds9 by storing the frame in shared memory, and thus ds9 updates
    without any additional delays and without any I/O overheads (see the :kbd:`--shmem` option below); this may be advantageous in cases where input and
    output files are stored on a network drive, or on a machine where file-caching is limited. Downside of this
    approach is that ds9 does not update the intensity histograms, and therefore automatic scaling (e.g. zscale) may
    not work (manually setting the intensity scaling still works fine though).


To run, one needs to specify a directory to monitor, a staging dir, and some
optional options that modify some of the internal algorithm::

    nirwals_watchdog incoming_data_dir --nonlinearity=/some/where/nonlinpoly.fits \
       --stage=/some/other/dir


Available options
==================

A full listing of all supported options can be obtained by running::

    nirwals_reduce --help

::

    ...> nirwals_watchdog --help
    usage: nirwals_watchdog [-h] [--stage STAGING_DIR] [--nonlinearity NONLINEARITY_FN]
                            [--refpixel REF_PIXEL_MODE] [--test TEST] [--nowait] [--shmem]
                            [--write] [--every EVERY] directory

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
      --every EVERY         only process every N-th frame instead of all frames


Options mirror the naming convention and functionality of ``nirwals_reduce`` wherever possible.


:kbd:`--stage`
  staging directory. This is where reduced files are written to before being sent to ds9 for
  displaying.

:kbd:`--nonlinearity`
  specify the non-linearity correction parameter file (see documentation for
  ``nirwals_reduce`` for more details)

:kbd:`--refpixel`
  reference pixels mode (use **blockyslope2** for best performance; see ``nirwals_reduce`` for
  further details)

:kbd:`--test`
  Run a test, simulating the arrival of new files in the monitored
  directory.

  Syntax is :kbd:`--test=delay:@filelist`, with delay giving a delay time
  between successive frames in seconds, and filelists specifying a file containing
  a list of "newly arrived" files, with one file per line (lines starting with #
  are ignored)

:kbd:`--nowait`
  By default, :command:`nirwals_watchdog` waits for a valid SAMP connection before starting work. If this
  is not required and/or desired, this waiting can be disabled. In that case, all on-the-fly reduction still works,
  but updating frames in ds9 is disabled.

  Note also that nirwals_watchdog only checks for a valid SAMP connection
  during startup, so loosing a SAMP connection during execution may result in an error, and starting up a SAMP
  server after the initial waiting phase will not automatically establish a connection.

:kbd:`--shmem`
  Load frames into ds9 using shared memory rather than writing files to disk and commanding ds9 to load
  the file from disk. See note above for strategies and implications.

:kbd:`--write`
  By default, reduced files are written to disk, unless the :kbd:`--shmem` option is selected, in which
  case writing files to disk is not necessary for proper watchdog operation. To force writing all reduced files
  to disk even if not strictly required, use this :kbd:`--write` flag.

:kbd:`--every=N`
  Instead of reducing every single frame as it comes in, limit processing to every N-th frame only. Note that this
  limits the frames being processed, not only displayed.


Run-times were tested on a modern laptop (i7 CPU, 32 GB RAM). Using a "full"
reduction mode, including non-linearity takes approx 0.4 to 0.5 seconds per frame,
from finding the newly arrived frame to end of writing the final result file.
Given the minimum read time of the NIRWALS instrument of ~0.7 seconds this should allow
monitoring incoming data in effectively real time (i.e. the previous frame is displayed
before the next read is fully read out).

In the case of using the `--every` command, the output is modified to look like this (here also using the --test
option)::

    ...> nirwals_watchdog raw/ --test=1:@N202310140001.list --stage=stage/ --every=3
    NirwalsWatchdog: Connecting to ds9 via SAMP protocol
    NirwalsWatchdog: Successfully connected to SAMPhub
    NirwalsWatchdog: Writing reduced files to forward to ds9
    NirwalsWatchdog: Starting on-the-fly reduction process
    NirwalsWatchdog: Starting to watch directory for new files: /nas/t7black/salt/SALT_data_2023/1014/raw
    NirwalsWatchdog: Test-mode activated, feeding 824 files with delays of 1.00 seconds
    WatchdogProcess: Starting new sequence: raw/N202310140001.1.1.1.fits
    WatchdogProcess: Sequence start is GOOD!
    WatchdogProcess: On-the-fly processing of raw/N202310140001.1.1.1.fits --> NEW SEQ completed after 0.061 seconds
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.2.fits), currently on frame 1 of 3
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.3.fits), currently on frame 2 of 3
    WatchdogProcess: On-the-fly processing of raw/N202310140001.1.1.4.fits --> N202310140001.1.1.4__qred.fits completed after 0.144 seconds
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.5.fits), currently on frame 1 of 3
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.6.fits), currently on frame 2 of 3
    WatchdogProcess: On-the-fly processing of raw/N202310140001.1.1.7.fits --> N202310140001.1.1.7__qred.fits completed after 0.220 seconds
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.8.fits), currently on frame 1 of 3
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.9.fits), currently on frame 2 of 3
    WatchdogProcess: On-the-fly processing of raw/N202310140001.1.1.10.fits --> N202310140001.1.1.10__qred.fits completed after 0.157 seconds
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.11.fits), currently on frame 1 of 3
    WatchdogProcess: Ignoring new frame (N202310140001.1.1.12.fits), currently on frame 2 of 3
    WatchdogProcess: On-the-fly processing of raw/N202310140001.1.1.13.fits --> N202310140001.1.1.13__qred.fits completed after 0.153 seconds

