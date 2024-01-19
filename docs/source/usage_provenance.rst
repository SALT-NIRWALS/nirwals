****************************************************************
Reading data provenance information using *nirwals_provenance*
****************************************************************

Each generated output file includes a provenance table containing detailed information about each input file,
all calibration products used during the processing, as well as information about the software environment in
which the file was generated. This information can be presented directly after the processing is done as part of
the execution of ``nirwals_reduce``, but can also be read back at any later time using the ``nirwals_provenance``
tool.

To read the information, run::

    ..> nirwals_provenance some_output_reduced.fits

For an example file, this call generates the following output::

    > nirwals_provenance N202401170002.2.2.reduced.fits
    DataProvenance:
     ==== DATA PROVENANCE INVENTORY ====
              invocation: /home/rss/.pyenv/versions/nirwals/bin/nirwals_reduce \
                            --nonlinearity=../nonlinpoly_N202310140001_1_1__new.fits --dumps=all \
                            --refpixel=blockyslope2 raw/N202401170002.2.2.1.fits
         python-version:: 3.12.1
         python-compiler: GCC 7.5.0
            python-build: ('main', 'Jan 18 2024 17:49:16')
              os-version: Linux-5.3.18-lp152.66-default-x86_64-with-glibc2.26
              os-version: #1 SMP Tue Mar 2 13:18:19 UTC 2021 (73933a3)
                os-uname: Linux stereo 5.3.18-lp152.66-default #1 SMP Tue Mar 2 13:18:19 \
                                   UTC 2021 (73933a3) x86_64 x86_64
               os-system: Linux
                 os-node: stereo
              os-release: 5.3.18-lp152.66-default
              os-machine: x86_64
            os-processor: x86_64
             interpreter: 64bit ELF
         nirwals-version: 0.1.0
         astropy-version: 6.0.0
           numpy-version: 1.26.3
           scipy-version: 1.11.4
         socket-hostname: stereo
                username: rss
              ref-header: /nas/zfspool_36x3/rss/SALT_data_2024/0117/raw/N202401170002.2.2.1.fits
                   input: /nas/zfspool_36x3/rss/SALT_data_2024/0117/raw/N202401170002.2.2.1.fits
                   input: /nas/zfspool_36x3/rss/SALT_data_2024/0117/raw/N202401170002.2.2.2.fits
    reference-pixel-mode: blockyslope2
           non-linearity: /nas/zfspool_36x3/rss/SALT_data_2024/nonlinpoly_N202310140001_1_1__new.fits
           urg-algorithm: linreg
     ==== DATA PROVENANCE INVENTORY END ====

This should allow to precisely regenerate the output file from its raw data even at a later time, thus ensuring
reproducibility of all generated data products.

For more details, you can also run ``nirwals_provenance --help`` to obtain the following usage information::

    ....> nirwals_provenance --help
    usage: nirwals_provenance [-h] files [files ...]

    positional arguments:
      files       list of input filenames

    options:
      -h, --help  show this help message and exit

