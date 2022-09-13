
import logging
import os

import astropy.io.fits as pyfits

class DataProvenance( object ):

    def __init__(self, logger=None, track_machine_data=False):

        if (logger is None):
            logger = logging.getLogger("DataProvenance")
        self.logger = logger

        self.inventory = {}

        if (track_machine_data):

            try:
                import platform
                self.add("python-version:", str(platform.python_version()))
                self.add("python-compiler", str(platform.python_compiler()))
                self.add("python-build", str(platform.python_build()))
                self.add("os-version", str(platform.platform()))

                self.add("os-uname", " ".join(platform.uname()))
                self.add("os-system", str(platform.system()))
                self.add("os-node", str(platform.node()))
                self.add("os-release", str(platform.release()))
                self.add("os-version", str(platform.version()))
                self.add("os-machine", str(platform.machine()))
                self.add("os-processor", str(platform.processor()))

                self.add("interpreter", " ".join(platform.architecture()))
            except:
                self.logger.debug("platform information not available for provenance")
                pass

            try:
                import socket
                self.add("socket-hostname", socket.gethostname())
            except:
                self.logger.debug("socket info not available, missing package socket")
                pass

            try:
                import getpass
                self.add("username", getpass.getuser())
            except:
                self.logger.debug("username not available, missing package getpass")
                pass


    def add(self, record, value):

        # check if this data-type is already given in the inventory
        if (record not in self.inventory):
            self.logger.debug("Adding new datatype to inventory: %s" % (record))
            self.inventory[record] = []

        self.logger.debug("Adding to data provenance inventory [%s]: %s" % (
            record, value
        ))

        if (os.path.isfile(value) or os.path.isdir(value)):
            values = os.path.abspath(value)

        self.inventory[record].append(value)

        return


    def write_to_header(self, hdr):
        return


    def write_as_hdu(self):
        return

    def report(self):
        print("\n ==== DATA PROVENANCE INVENTORY ==== ")
        for key in self.inventory:
            for fn in self.inventory[key]:
                print("% 20s: %s" % (key, fn))
        print(" ==== DATA PROVENANCE INVENTORY END ==== ")
        return

