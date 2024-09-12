import logging

import numpy


class NirwalsGainData( object ):
    logger = logging.getLogger("NirwalsGainData")
    gains = numpy.array([2.0] * 32)
    name = "Fixed_Gain2_Fallback"

    def __init__( self):
        self.logger.debug("Setting up NirwalsGainData (%s) : %s" % (
            self.name, ",".join(["%.3f" % g for g in self.gains])
            )
        )

class NirwalsGainData_GainOne( NirwalsGainData ):
    gains = numpy.array([1.0] * 32)
    name = "No_Gain_Correction"

class NirwalsGainData_Mean( NirwalsGainData ):
    name  = "NirwalsGainData_Mean"

    def __init__( self, gain_data ):
        self.logger.debug("Generating MEAN gain values from %s" % (gain_data.name))
        _gains = gain_data.gains
        average_gain = numpy.nanmean(_gains)
        self.gains = numpy.full_like(_gains, average_gain)
        self.name = gain_data.name+"__Mean"
        super(NirwalsGainData_Mean, self).__init__()


class NirwalsGainData__Commissioning( NirwalsGainData ):
    name = "Commissioning"
    # To check: Amp 18,22,26 (16, 4?)
    gains = numpy.array([2.00505, 1.91840, 1.92488, 1.94229, 1.87510, 1.97130, 1.96191, 1.93191, #  0- 7
                         1.86682, 1.82613, 1.91389, 2.10541, 2.03704, 2.05165, 2.02557, 1.98288, #  8-15
                         1.78187, 1.93506, 1.53476, 1.88753, 1.94853, 1.86109, 1.64050, 1.99922, # 16-23
                         1.60471, 1.97976, 1.69511, 1.91282, 1.87358, 1.91070, 1.96332, 1.93810, # 24-32
                         ])


class NirwalsGainData__ByDate( NirwalsGainData) :

    def __init__(self, date_obs):
        self.logger.debug("Checking NirwalsGainData by Date (%s)" % (date_obs))
        # add some more code here once more values are found
        _gain = NirwalsGainData__Commissioning()

        self.name = "ByDate__" + _gain.name
        self.gains = _gain.gains
        super(NirwalsGainData__ByDate, self).__init__()






class NirwalsGain( object ):

    gain_data = None
    date_obs = None
    amp_width = 64 # each gain is 64 pixels wide
    logger = logging.getLogger("NirwalsGain")

    def __init__(self, date_obs=None, header=None, gain_mode='none'):

        self.gain_mode = gain_mode

        if (header is not None):
            # if we have a valid FITS header we can look up the date from there
            try:
                self.date_obs = header['DATE-OBS']
            except:
                self.logger.warning("Found Header, but no DATE-OBS")
                pass

        if (date_obs is not None):
            self.date_obs = date_obs

        self.logger.info("USING GAIN MODE: %s (date: %s)" % (self.gain_mode, self.date_obs))


        # No gain correction requested
        if (self.gain_mode == 'none'):
            self.gain_data = NirwalsGainData_GainOne()

        # We don't know what date, so we'll use a sensible default
        elif (self.date_obs is None):
            self.gain_data = NirwalsGainData()

        # If plain correction --> Use the default value
        elif (self.gain_mode == 'plain'):
            self.gain_data = NirwalsGainData()

        # full correction: Determine what set of GAIN values to use, based on DATE-OBS
        elif (self.gain_mode == 'full' and self.date_obs is not None):
            self.gain_data = NirwalsGainData__ByDate(self.date_obs)

        # averaged gain values --> look up detailed values, then average
        elif (self.gain_mode == 'average'):
            if (self.date_obs is None):
                _gain_data = NirwalsGainData()
            else:
                _gain_data = NirwalsGainData__ByDate(self.date_obs)
            self.gain_data = NirwalsGainData_Mean(_gain_data)

        # in all other cases: also use the default values
        else:
            self.gain_data = NirwalsGainData()


    def get_name(self):
        if (self.gain_data is not None):
            return self.gain_data.name
        return "UNKNOWN"

    def by_amp(self, amp):
        try:
            gain = self.gain_data.gains[amp]
        except IndexError:
            gain = 2.0
        return gain

    def get_gains(self):
        return self.gain_data.gains

    def amp_corrections(self):
        full = self.gain_data.gains.repeat(self.amp_width)
        return full



