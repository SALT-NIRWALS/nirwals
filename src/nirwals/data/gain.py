
import numpy


class NirwalsGainData( object ):
    gains = numpy.array([2.0] * 32)
    name = "Fixed_Gain2_Fallback"


class NirwalsGainData__Commissioning( NirwalsGainData ):
    name = "Commissioning"
    # To check: Amp 18,22,26 (16, 4?)
    gains = numpy.array([2.00505, 1.91840, 1.92488, 1.94229, 1.87510, 1.97130, 1.96191, 1.93191, #  0- 7
                         1.86682, 1.82613, 1.91389, 2.10541, 2.03704, 2.05165, 2.02557, 1.98288, #  8-15
                         1.78187, 1.93506, 1.53476, 1.88753, 1.94853, 1.86109, 1.64050, 1.99922, # 16-23
                         1.60471, 1.97976, 1.69511, 1.91282, 1.87358, 1.91070, 1.96332, 1.93810, # 24-32
                         ])
    pass


class NirwalsGain( object ):

    gain_data = None
    date_obs = 0
    amp_width = 64 # each gain is 64 pixels wide

    def __init__(self, date_obs=None, header=None):

        if (header is not None):
            # if we have a valid FITS header we can look up the date from there
            try:
                self.date_obs = header['DATE-OBS']
            except:
                pass

        if (date_obs is None):
            self.date_obs = date_obs

        if (self.date_obs is None):
            # we don't know anything about this frame, just use some default values
            self.gain_data = NirwalsGainData()
        else:
            # Determine what set of GAIN values to use, based on DATE-OBS
            self.gain_data = NirwalsGainData__Commissioning()


    def get_name(self):
        return self.gain_data.name

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



