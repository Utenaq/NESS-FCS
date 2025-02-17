# simulation molecular trajectories and single-molecule fluorescence trace

from trajectory import trajectory
from fluorescence_surf import fluorescence_surf
from reaction_3state import reaction_3state
#import minpy.numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import axes3d
import time
import sys
import os

# paramaters initials
dt = 10  # time interval (us)

# simulation time (us)
# 2**21: ~2 s
# 2**20: ~1 s
# 2**19: ~0.5 s
# 2**18: ~0.25 s
# 2**17: ~0.125 s
totalTime = 1000000

repeatCycle = 100  # repeat totalTime simulation
fileCycle = 1  # repeat simulation files

# simulatoin condition initialization
border = [3, 3, 15]
# moleculeNum =  np.int(sys.argv[3])*1000/12 #simulated molecule number
surfMoleculeNum = 10

pA = 1
pB = 0

AmoleculeNum = int(surfMoleculeNum * pA)
BmoleculeNum = int(surfMoleculeNum * pB)
CmoleculeNum = surfMoleculeNum - AmoleculeNum - BmoleculeNum


# fluor quantum yield
Qfluor = 1
QA = np.float(sys.argv[7])
QB = np.float(sys.argv[8])
QC = np.float(sys.argv[9])
QacceptorA = 0
QdonorB = 0
QacceptorB = 0

# reaction detail
# Input: k01 k10 k12 k21 k20 k02
k_Matrix = np.zeros((3, 3))
k_Matrix[0][1] = np.int(sys.argv[1]) / 1e6  # us-1
k_Matrix[1][0] = np.int(sys.argv[2]) / 1e6  # us-1
k_Matrix[1][2] = np.int(sys.argv[3]) / 1e6  # us-1
k_Matrix[2][1] = np.int(sys.argv[4]) / 1e6  # us-1
k_Matrix[2][0] = np.int(sys.argv[5]) / 1e6  # us-1
k_Matrix[0][2] = np.int(sys.argv[6]) / 1e6  # us-1
k_Matrix[0][0] = -k_Matrix[0][1] - k_Matrix[0][2]
k_Matrix[1][1] = -k_Matrix[1][0] - k_Matrix[1][2]
k_Matrix[2][2] = -k_Matrix[2][0] - k_Matrix[2][1]

initInfo = 'trajectory simulation time interval: ' + str(dt) + ' us\n' + \
           'total simulation time: ' + str(totalTime * 1e-6) + 'x' + str(repeatCycle) + ' s\n' + \
           'molecule number: ' + str(surfMoleculeNum) + '\n' + \
           'border volume: ' + 'x'.join(map(str, border)) + ' um^3\n' + \
           'apparent concentration: ' + str(
    int(surfMoleculeNum / 6.023 / border[0] / border[1] / border[2] * 1e4)) + ' pM\n' + \
           '-------------------------------------\n' + \
           'A ratio: ' + str(pA) + '\n' + \
           'B ratio: ' + str(pB) + '\n' + \
           'A molecule number: ' + str(AmoleculeNum) + '\n' + \
           'B molecule number: ' + str(BmoleculeNum) + '\n' + \
           'C molecule number: ' + str(CmoleculeNum) + '\n' + \
           'A effective brightness: ' + str(Qfluor * 4.172) + '\n' + \
           '-------------------------------------\n\n' + \
           'Reaction rate constant : ' + (' ').join(map(lambda x: (' ').join(map(str, x)), k_Matrix)) + '\n' + \
           'Intensity A:' + str(QA) + '\n' + \
           'Intensity B:' + str(QB) + '\n' + \
            'Intensity C:' + str(QC) + '\n' + \
           '-------------------------------------\n\n' + \
           '1-typical, 2-without poisson random, 3-typical without int, 4-without poisson random & int' + \
           'without poisson noise & crosstalk.\n\n'

# log information file;

path = '../LnrSimulation/180524serie/' + 'K_'+('_').join(map(lambda x: ('_').join(map(str, x)), k_Matrix)) + '_QA' + str(QA) + '_QB' + str(QB) + '_QC' + str(QC) + '_pA' + str(
    pA) + '_pB' + str(pB)

try: os.makedirs(path)
except: pass

with open(path + '/log.txt', 'w') as f:
    f.write(initInfo + str(time.asctime()) + ' - trace: simulating...\n\n')

# plot
#plt.ion()
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
x = np.arange(0, totalTime*1e-3, dt*1e-3)
ax1.set_xlabel('t/ ms')
ax2.set_xlabel('t/ ms')
line1, = ax1.plot(x,np.zeros(np.int(totalTime/dt)), c='r')
line2, = ax1.plot(x,np.zeros(np.int(totalTime/dt)), c='b')
line31, = ax2.plot(x,np.zeros(np.int(totalTime/dt)))
line32, = ax2.plot(x,np.zeros(np.int(totalTime/dt)))
line33, = ax2.plot(x,np.zeros(np.int(totalTime/dt)))
'''
fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
x_m = np.arange(0, 1000*dt*1e-3, dt*1e-3)
ax1.set_xlabel('t/ ms')
ax2.set_xlabel('t/ ms')
line101, = ax1.plot(x_m,np.zeros(1000), c='r')
line102, = ax1.plot(x_m,np.zeros(1000), c='b')
line211, = ax2.plot(x_m,np.zeros(1000))
line212, = ax2.plot(x_m,np.zeros(1000))
line213, = ax2.plot(x_m,np.zeros(1000))
'''

# ax2.plot(x, reactionTrajectory.moleculenum[1][:-10000], c='g')
# ax2.plot(x, reactionTrajectory.moleculenum[2][:-10000], c='b')

def update(t,data01,data02,data11):
    ax1.set_ylim(0,30)
    ax2.set_ylim(0,10)
    line1.set_data (t, data01)
    line2.set_data (t, data02)
    line31.set_data (t, data11[:, 0])
    line32.set_data(t, data11[:, 1])
    line33.set_data(t, data11[:, 2])
    #plt.pause(0.001)
    #ax1.figure.canvas.draw()
    #ax2.figure.canvas.draw()
def update_monitor(t,data01,data02,data11):
    ax1.set_ylim(0,20)
    ax2.set_ylim(0,20)
    line1.set_data (t, data01)
    line2.set_data (t, data02)
    line31.set_data (t, data11[:, 0])
    line32.set_data(t, data11[:, 1])
    line33.set_data(t, data11[:, 2])

for fileNum in range(fileCycle):
    # initial stochastic coordinates generation
    initPosition = [np.random.uniform(-border[0] / 2, border[0] / 2, surfMoleculeNum),
                    np.random.uniform(-border[1] / 2, border[1] / 2, surfMoleculeNum),
                    np.repeat(-border[2] / 2, surfMoleculeNum)]
    initState = np.zeros(surfMoleculeNum)
    for i in range(surfMoleculeNum): initState[i] = QA if i < AmoleculeNum else (
        QB if (i >= AmoleculeNum and i < AmoleculeNum + BmoleculeNum) else QC)

    # donor channel trace file;
    # accptor channel trace file;
    with open(path + '/donor_' + str(fileNum) + '.dat', 'a') as f:
        f.write('')

    with open(path + '/donor_' + str(fileNum) + 'nr.dat', 'a') as f:
        f.write('')

    with open(path + '/moleculenum_' + str(fileNum) + '.dat', 'a') as f:
        f.write('')

    for repeatNum in range(repeatCycle):
        print('Simulating: Cycle' + str(repeatNum) + "\n")
        # trajectory simulation

        # reaction simulation
        reactionTrajectory = reaction_3state(dt, totalTime, surfMoleculeNum, initState)
        reactionTrajectory.react(k_Matrix, QA, QB, QC)
        initState = reactionTrajectory.state[-1, :]
        fluoreDonor = fluorescence_surf(dt, Qfluor)
        fluoreDonor.collectPhoton(reactionTrajectory.state)

        # collection fluorescence
        # Donor channel
        ##Acceptor channel
        # fluoreAcceptor = fluorescence(dt, Qfluor = Qacceptor)
        # fluoreAcceptor.collectPhoton(molecularTrajectory.positionX, molecularTrajectory.positionY, molecularTrajectory.positionZ)
        #print(np.shape(reactionTrajectory.moleculeNum))
        update(x,fluoreDonor.trace,fluoreDonor.trace_nr,reactionTrajectory.moleculenumtrace)
        if repeatNum==0: plt.savefig(path+'/InitTrace.png')
        #ax2.plot(x, reactionTrajectory.moleculenumtrace[:,0])
        #plt.ioff()
        #plt.show()
        with open(path + '/donor_' + str(fileNum) + '.dat', 'a') as f:
            np.savetxt(f, fluoreDonor.trace, fmt='%.3f')

        with open(path + '/donor_' + str(fileNum) + 'nr.dat', 'a') as f:
            np.savetxt(f, fluoreDonor.trace_nr, fmt='%.3f')

        with open(path + '/moleculenum_' + str(fileNum) + '.dat', 'a') as f:
            np.savetxt(f, reactionTrajectory.moleculenumtrace, fmt='%.3f')
            # np.savetxt(f, fluoreDonor.singleTrace, fmt='%i')
        # with open('./acceptor_'+str(fileNum)+'.txt', 'a') as f:
        #    np.savetxt(f, fluoreAcceptor.trace, fmt='%i')
        #    #np.savetxt(f, fluoreAcceptor.singleTrace, fmt='%i')
#plt.ioff()
plt.savefig(path+'/EndTrace.png')
#plt.show()
with open(path + '/log.txt', 'a') as f:
    f.write(str(time.asctime()) + ' - trace: done.\n\n')
