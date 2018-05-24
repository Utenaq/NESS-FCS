# trajectories of molecules gereation
# Wiener process (Brownian motion)
# peroidic boundary condition

#import minpy.numpy as np
import numpy as np

class reaction_3state:
    """gerenate molecule trajectories"""

    def __init__(self, dt, totalTime, moleculeNum, initState):
        self.dt = dt  # time intervel (us)
        self.totalTime = totalTime  # length of simulated time (us)
        self.moleculeNum = moleculeNum  # number of simulated molecules
        self.stepNum = np.int(self.totalTime/self.dt)  # simulation steps
        self.size = self.stepNum + 1000000
        self.state3 = np.zeros([self.size, self.moleculeNum])
        self.moleculenum3 = np.zeros([self.size, 3])
        self.state3[0, :] = initState

    def react(self, k_Matrix, QA, QB, QC):
        # state=i, kplus:1->0, kminus:0->1
        for i in np.arange(self.moleculeNum):
            j = 0
            while (j < self.stepNum):
                if self.state3[j, i] == QA: #state 0
                    time = np.int(np.round(np.random.exponential(1 / (k_Matrix[0][1]+k_Matrix[0][2])) / self.dt))
                    if time > 0:
                        self.state3[j:j + time, i] = self.state3[j, i]
                        self.moleculenum3[j:j + time, 0] += 1
                        j = j + time
                    jump_judge = k_Matrix[0][1] / (k_Matrix[0][1]+k_Matrix[0][2])
                    if np.random.uniform(0,1)<jump_judge:
                        self.state3[j, i] = QB
                    else:
                        self.state3[j, i] = QC
                if self.state3[j, i] == QB: #state 1
                    time = np.int(np.round(np.random.exponential(1 / (k_Matrix[1][2]+k_Matrix[1][0])) / self.dt))
                    if time > 0:
                        self.state3[j:j + time, i] = self.state3[j, i]
                        self.moleculenum3[j:j + time, 1] += 1
                        j = j + time
                    jump_judge = k_Matrix[1][2] / (k_Matrix[1][2]+k_Matrix[1][0])
                    if np.random.uniform(0,1)<jump_judge:
                        self.state3[j, i] = QC
                    else:
                        self.state3[j, i] = QA
                if self.state3[j, i] == QC: #state 2
                    time = np.int(np.round(np.random.exponential(1 / (k_Matrix[2][0]+k_Matrix[2][1])) / self.dt))
                    if time > 0:
                        self.state3[j:j + time, i] = self.state3[j, i]
                        self.moleculenum3[j:j + time, 2] += 1
                        j = j + time
                    jump_judge = k_Matrix[2][0] / (k_Matrix[2][0]+k_Matrix[2][1])
                    if np.random.uniform(0,1)<jump_judge:
                        self.state3[j, i] = QA
                    else:
                        self.state3[j, i] = QB
            self.state = self.state3[0:self.stepNum, :]
            self.moleculenumtrace = self.moleculenum3[0:self.stepNum, :]