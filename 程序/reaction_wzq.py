#trajectories of molecules gereation
#Wiener process (Brownian motion)
#peroidic boundary condition
    
import numpy as np

class reaction_wzq:
    """gerenate molecule trajectories"""
    def __init__(self, dt, totalTime, moleculeNum):
        self.dt = dt    #time intervel (us)
        self.totalTime = totalTime    #length of simulated time (us)
        self.moleculeNum = moleculeNum    #number of simulated molecules
        self.stepNum = self.totalTime/self.dt    #simulation steps
        self.size = self.stepNum + 10000
        self.state2 = np.zeros([self.size, self.moleculeNum])


    def react(self, kplus, kminus, Q):
        #bright==1 dark==0, kplus:1->0, kminus:0->1
        initstate_judge = kplus/(kplus+kminus)
        for i in np.arange(self.moleculeNum):
            if np.random.uniform(0,1)<initstate_judge:
                self.state2[0,i] = Q
            else:
                self.state2[0,i] = 1
            j = 0
            while (j<self.stepNum):
                if self.state2[j,i]==1:
                    time = np.round(np.random.exponential(1/kplus)/self.dt)
                    if time>0:
                        self.state2[j:j+time,i] = self.state2[j,i]
                        j = j+time
                    self.state2[j,i] = Q
                if self.state2[j,i]==Q:
                    time = np.round(np.random.exponential(1/kminus)/self.dt)
                    if time>0:
                        self.state2[j:j+time,i] = self.state2[j,i]
                        j = j+time
                    self.state2[j,i] = 1
            self.state = self.state2[0:self.stepNum,:]