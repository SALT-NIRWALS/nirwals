************
Installation
************

.. _install:

Installing ``nirwals``
**********************


Using pip
=========

.. Warning::
   Since nirwals is still actively being developed, we recommend to follow the
   instructions for installing the latest version of nirwals :ref:'directly
   from github <install_from_github>'


.. _install_from_github:
Get latest version from Github
================================

To download the latest version from github::

        git clone  https://github.com/SALT-NIRWALS/nirwals.git
        # that should put all software into a nirwals directory

Once completed, install using pip::

        cd nirwals
        pip install .

In addition to the nirwals package itself, this will install all the required packages.



Requirements
***************

``nirwals`` has the following requirements:

- `Python`, minimum version is 3.8

- `numpy`

- `scipy`

- `matplotlib`

- `pandas`

- `multiparlog`

For the watchdog to work, the following packages are required:

- `pyvo`

- `watchdog`

- `sysv_ipc` (this is strictly only required for faster display of reduced data in ds9)
