#fluorescence simulation

import numpy as np

class fluorescence_wzq:
    """fluorescence emision simulation"""
    def __init__(self, dt, Qfluor, power = 100, wavelength = 532, radius = 0.3, axialRatio = 5, Qdetect = 0.01, sigma = 2.2e-8):
        self.dt = dt #time interval (s)
        self.power = power #laser power (uW)
        ePhoton = {
                '488':4.07e-19,
                '532':3.73e-19,
                '633':3.14e-19
                }
        self.ePhoton = ePhoton[str(wavelength)] #energy of single photon (J)
        self.radius = radius #distance from axial where laser intensity drops by 1/e^2 (um)
        self.axialRatio = axialRatio
        self.Qdetect = Qdetect #detection efficiency
        self.sigma = sigma #absorption cross section (um^2)
        self.Qfluor = Qfluor

    def collectPhoton(self, x, y, z, react):

        self.x = x
        self.y = y
        self.z = z/self.axialRatio #rescale of z coordinate
        self.react=react

        self.stepNum = len(x)
        self.moleculeNum = x.size/len(x)
        
        #print 'calculating fluorescence...'
        self.Imax = 2*self.power/(np.pi*self.radius**2) # the center laser intensity (uW/um^2)
        self.photonRate = self.Imax*self.Qfluor*self.react/self.ePhoton*self.Qdetect*self.sigma*1e-12 #(us^-1)
        
        self.singleTrace = np.random.poisson(self.photonRate*self.dt)*np.exp(-2*(self.x**2+self.y**2+self.z**2)/self.radius**2)
        self.trace = np.sum(self.singleTrace, axis=1)
        
        self.singleTrace_nr = self.photonRate*self.dt*np.exp(-2*(self.x**2+self.y**2+self.z**2)/self.radius**2)
        self.trace_nr = np.sum(self.singleTrace_nr, axis=1)
        #print 'suming all molecular fluorescence at the same time'
               
        
        
