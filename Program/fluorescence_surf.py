# fluorescence simulation

#import minpy.numpy as np
import numpy as np


class fluorescence_surf:
    """fluorescence emission simulation"""

    def __init__(self, dt, Qfluor, power=100, wavelength=532, radius=0.3, axialRatio=5, Qdetect=0.01, sigma=2.2e-8):
        self.dt = dt  # time interval (s)
        self.power = power  # laser power (uW)
        ePhoton = {
            '488': 4.07e-19,
            '532': 3.73e-19,
            '633': 3.14e-19
        }
        self.ePhoton = ePhoton[str(wavelength)]  # energy of single photon (J)
        self.radius = radius  # distance from axial where laser intensity drops by 1/e^2 (um)
        self.axialRatio = axialRatio
        self.Qdetect = Qdetect  # detection efficiency
        self.sigma = sigma  # absorption cross section (um^2)
        self.Qfluor = Qfluor

    def collectPhoton(self, react):
        self.react = react

        # print 'calculating fluorescence...'
        self.I = self.power  # the center laser intensity (uW/um^2)
        self.photonRate = self.I * self.Qfluor * self.react / self.ePhoton * self.Qdetect * self.sigma * 1e-12  # (us^-1)

        self.singleTrace = np.random.poisson(self.photonRate * self.dt)
        self.trace = np.sum(self.singleTrace, axis=1)

        self.singleTrace_nr = self.photonRate * self.dt
        self.trace_nr = np.sum(self.singleTrace_nr, axis=1)
        # print 'suming all molecular fluorescence at the same time'



