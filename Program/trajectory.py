#trajectories of molecules gereation
#Wiener process (Brownian motion)
#peroidic boundary condition
    
#import minpy.numpy as np
import numpy as np

class trajectory:
    """gerenate molecule trajectories"""
    def __init__(self, dt, totalTime, moleculeNum, initPosition, border):
        self.dt = dt    #time intervel (us)
        self.totalTime = totalTime    #length of simulated time (us)
        self.moleculeNum = moleculeNum    #number of simulated molecules
        self.stepNum = np.int(self.totalTime/self.dt)    #simulation steps
        self.positionX = np.zeros([self.stepNum, self.moleculeNum])
        self.positionY = np.zeros([self.stepNum, self.moleculeNum])
        self.positionZ = np.zeros([self.stepNum, self.moleculeNum])
        self.border = np.array(border) #simulation box border length (um)
        #print(self.border)
        self.initPositionX = initPosition[0]
        self.initPositionY = initPosition[1]
        self.initPositionZ = initPosition[2]

    def diffusion(self, diffCoef, drift = 0):
        #diffCoef: diffusion coefficient (um^2/us), usually in the magnitude of 1e-4
        walkDist = np.sqrt(2*diffCoef*self.dt)
        self.positionX[0, :] = np.random.normal(np.repeat(drift*self.dt, self.moleculeNum), walkDist)+self.initPositionX
        self.positionY[0, :] = np.random.normal(np.repeat(drift*self.dt, self.moleculeNum), walkDist)+self.initPositionY
        self.positionZ[0, :] = np.random.normal(np.repeat(drift*self.dt, self.moleculeNum), walkDist)+self.initPositionZ
        for step in np.arange(self.stepNum-1)+1:
            self.positionX[step, :] = np.random.normal(np.repeat(drift*self.dt, self.moleculeNum), walkDist)+self.positionX[step-1, :]
            self.positionY[step, :] = np.random.normal(np.repeat(drift*self.dt, self.moleculeNum), walkDist)+self.positionY[step-1, :]
            self.positionZ[step, :] = np.random.normal(np.repeat(drift*self.dt, self.moleculeNum), walkDist)+self.positionZ[step-1, :]
            #print 'walking '+str(step*100/self.stepNum)+'% steps...'
            #periodic boundary conditoin         
            indX = np.abs(self.positionX[step, :]) > self.border[0]/2
            indY = np.abs(self.positionY[step, :]) > self.border[1]/2
            indZ = np.abs(self.positionZ[step, :]) > self.border[2]/2
            #symmetrical coordinate
            #x' = x +/- border_length
            if np.sum(indX):
                self.positionX[step, indX] = self.positionX[step, indX]-self.positionX[step, indX]/np.abs(self.positionX[step, indX])*self.border[0]
            if np.sum(indY):
                self.positionY[step, indY] = self.positionY[step, indY]-self.positionY[step, indY]/np.abs(self.positionY[step, indY])*self.border[1]
            if np.sum(indZ):
                self.positionZ[step, indZ] = self.positionZ[step, indZ]-self.positionZ[step, indZ]/np.abs(self.positionZ[step, indZ])*self.border[2]
