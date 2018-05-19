# simulation molecular trajectories and single-molecule fluorescence trace

from Program.Bimolecular_Reaction import state, molecule, surfMolecule, volumeMolecule
#import minpy.numpy as np
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import axes3d

import time
import sys

# paramaters initials
dt = 1000  # time interval (us)

Avogadro_Constant=6.022140857e23
# simulation time (us)
# 2**21: ~2 s
# 2**20: ~1 s
# 2**19: ~0.5 s
# 2**18: ~0.25 s
# 2**17: ~0.125 s
totalTime = 1000000

repeatCycle = 10  # repeat totalTime simulation
fileCycle = 10  # repeat simulation files

# simulatoin condition initialization
border = np.array([20, 5, 1])
drift = np.array([1e-5 ,0,0])

concentration=1.66113 #nM/l
Ndensity=concentration*Avogadro_Constant*1e-24 #atomcounts/pl
dt_flowin=(1/Ndensity)/(border[1]*border[2]*drift[0])
# moleculeNum =  np.int(sys.argv[3])*1000/12 #simulated molecule number
volMoleculeNum = int(Ndensity*(border[0]*border[1]*border[2]))
surfMoleculeNum = 40

pA = np.float(sys.argv[11])
pB = np.float(sys.argv[12])

AmoleculeNum = int(volMoleculeNum * pA)
BmoleculeNum = int(volMoleculeNum * pB)
CmoleculeNum = volMoleculeNum - AmoleculeNum - BmoleculeNum

Adiffcoef = 2e-5 * 1125 / np.int(sys.argv[7])  # component A diffusion coefficient (um^2/us), ~10ms
Bdiffcoef = 2e-5 * 1125 / np.int(sys.argv[7])  # component B diffusion coefficient (um^2/us)
Cdiffcoef = 2e-5 * 1125 / np.int(sys.argv[7])  # component C diffusion coefficient (um^2/us)

# fluor quantum yield
Qfluor = 1
Q0 = np.float(sys.argv[8])
Q1 = np.float(sys.argv[9])
Q2 = np.float(sys.argv[10])
QacceptorA = 0
QdonorB = 0
QacceptorB = 0
rigidRadius=5e-3 #纳米量级
bindingRadius=1e-2
unbindingRadius=2.5e-2
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
           'molecule number: ' + str(volMoleculeNum) + '\n' + \
           'border volume: ' + 'x'.join(map(str, border)) + ' um^3\n' + \
           'apparent concentration: ' + str(
    int(volMoleculeNum / 6.023 / border[0] / border[1] / border[2] * 1e4)) + ' pM\n' + \
           '-------------------------------------\n' + \
           'A ratio: ' + str(pA) + '\n' + \
           'B ratio: ' + str(pB) + '\n' + \
           'A molecule number: ' + str(AmoleculeNum) + '\n' + \
           'B molecule number: ' + str(BmoleculeNum) + '\n' + \
           'C molecule number: ' + str(CmoleculeNum) + '\n' + \
           'A diffusion coefficient: ' + str(Adiffcoef * 1e-6) + ' m^2/s\n' + \
           'B diffusion coefficient: ' + str(Bdiffcoef * 1e-6) + ' m^2/s\n' + \
           'C diffusion coefficient: ' + str(Cdiffcoef * 1e-6) + ' m^2/s\n' + \
           'A effective brightness: ' + str(Qfluor * 4.172) + '\n' + \
           '-------------------------------------\n\n' + \
           'Reaction rate constant : ' + (' ').join(map(lambda x: (' ').join(map(str, x)), k_Matrix)) + '\n' + \
           'Intensity A:' + str(Q0) + '\n' + \
           'Intensity B:' + str(Q1) + '\n' + \
           '-------------------------------------\n\n' + \
           '1-typical, 2-without poisson random, 3-typical without int, 4-without poisson random & int' + \
           'without poisson noise & crosstalk.\n\n'

# log information file;

E=state(0, Q0, 0)
S_ES=state(1, Q1, 2e-5 * 1125 / np.int(sys.argv[7]))
P_EP=state(2, Q2, 2e-5 * 1125 / np.int(sys.argv[7]))
E.chemlink(S_ES)
E.chemlink(P_EP)
S_ES.chemlink(E)
S_ES.chemlink(P_EP)
P_EP.chemlink(E)
P_EP.chemlink(S_ES)
states=np.array([E, S_ES, P_EP])

path = '../NonLnrSimulation/180517serie/D' + sys.argv[7] + '_QA' + str(Q0) + '_QB' + str(Q1) + '_QC' + str(Q2) + '_pA' + str(
    pA) + '_pB' + str(pB)

with open(path + '/log.txt', 'w') as f:
    f.write(initInfo + str(time.asctime()) + ' - trace: simulating...\n\n')

# initial stochastic coordinates generation for volmol
    volInitPosition = np.array([np.random.uniform(-border[0] / 2, border[0] / 2, volMoleculeNum),
                                np.random.uniform(-border[1] / 2, border[1] / 2, volMoleculeNum),
                                np.random.uniform(-border[2] / 2, border[2] / 2, volMoleculeNum)]).transpose()

    volPosition=volInitPosition
    # initial state generation for volmol
    surfinitState = np.repeat(states[0], surfMoleculeNum)
    # initial stochastic coordinates generation for surfmol
    surfInitPosition = np.array([np.random.uniform(-border[0] / 2 + 2, border[0] / 2 - 2, surfMoleculeNum),
                                             np.random.uniform(-border[1] / 2, border[1] / 2, surfMoleculeNum),
                                             np.repeat(-border[2] / 2, surfMoleculeNum)]).transpose()
    surfPosition=surfInitPosition
    # initial state generation for surfmol
    volInitStateID= np.random.binomial(1, pB, size=volMoleculeNum) + np.repeat(1,volMoleculeNum)
    volInitState = states[volInitStateID]

    #generate Arrays of surfmol and volmol
    surfIDList=np.arange(1,surfMoleculeNum+1)
    volIDList = np.arange(1, volMoleculeNum + 1)

    allVolMolecules = [volumeMolecule(volIDList[j], volInitPosition[j], volInitState[j],None) for j in range (volMoleculeNum)]
    allSurfMolecules = [surfMolecule(surfIDList[i], surfInitPosition[i], surfinitState[i], None, rigidRadius,
        bindingRadius, unbindingRadius, k_Matrix) for i in range(surfMoleculeNum)]

    molID=volMoleculeNum

for fileNum in range(fileCycle):
    EnzymePositionArray=[]
    SPositionArray=[]
    PPositionArray=[]

    """
    surfPositionForPlot =surfPosition.transpose()  # position data
    x, y, z = surfPositionForPlot[0], surfPositionForPlot[1], surfPositionForPlot[2]
    surfPlot = plt.subplot(211)  # 创建一个三维的绘图工程
    #  将数据点分成三部分画，在颜色上有区分度
    for point in range(surfMoleculeNum):
        if allSurfMolecules[point].state == states[1]:
            surfPlot.scatter(x[point], y[point], c='r')
        if allSurfMolecules[point].state == states[2]:
            surfPlot.scatter(x[point], y[point], c='g')
        if allSurfMolecules[point].state == states[0]:
            surfPlot.scatter(x[point], y[point], c='b')

    surfPlot.set_ylabel('Y')
    surfPlot.set_xlabel('X')

    volPositionForPlot = volPosition.transpose()  # position data
    x, y, z = volPositionForPlot[0], volPositionForPlot[1], volPositionForPlot[2]
    volPlot = plt.subplot(212, projection='3d')  # 创建一个三维的绘图工程
    for point in range(volMoleculeNum):
        if allVolMolecules[point].state == states[1]:
            volPlot.scatter(x[point], y[point], z[point], c='r')
        if allVolMolecules[point].state == states[2]:
            volPlot.scatter(x[point], y[point], z[point], c='g')

    volPlot.set_zlabel('Z')  # 坐标轴
    volPlot.set_ylabel('Y')
    volPlot.set_xlabel('X')
    plt.show()
    """

    # donor channel trace file;
    # accptor channel trace file;
    with open(path + '/donor_' + str(fileNum) + '.txt', 'a') as f:
        f.write('')

    with open(path + '/donor_' + str(fileNum) + 'nr.txt', 'a') as f:
        f.write('')

    for repeatNum in range(repeatCycle):
        countsArray = np.zeros((np.int(totalTime / dt), 6))
        print('Simulating: Cycle' + str(repeatNum) + "\n")
        t_cycle=0
        t_flowin=dt_flowin
        while totalTime>t_cycle:
            #if flow
            # reaction-diffusion simulation
            # temporary position
            countsArray[np.int(t_cycle/dt)][0]=t_cycle
            tempEnzymePosition = [[],[],[]]
            tempSPosition = [[],[],[]]
            tempPPosition = [[],[],[]]
            # update E binding
            for a in allVolMolecules:
                for b in allSurfMolecules:
                    b.bireactTest(a,dt)
                #allSurfMolecules.bireactTest(allVolMolecules,dt)

            # update ES EP reaction
            # update transition time
            for b in allSurfMolecules:
                b.update(dt)
                if b.state == E:
                    countsArray[np.int(t_cycle / dt)][1] += 1
                elif b.state == S_ES:
                    countsArray[np.int(t_cycle / dt)][2] += 1
                elif b.state == P_EP:
                    countsArray[np.int(t_cycle / dt)][3] += 1
                for i in range(3):
                    tempEnzymePosition[i].append(b.getposition()[i])

            # diffusion
            for a in allVolMolecules:
                a.update(dt,border,drift)
                #print(a.state.ID==1)
                if a.state==S_ES:
                    for i in range (3): tempSPosition[i].append(a.getposition()[i])
                    countsArray[np.int(t_cycle / dt)][4] += 1
                elif a.state==P_EP:
                    for i in range(3): tempPPosition[i].append(a.getposition()[i])
                    countsArray[np.int(t_cycle / dt)][5] += 1
                if a.isinbox()==False:
                    allVolMolecules.remove(a)
            #data analysis
            EnzymePositionArray.append(np.array(tempEnzymePosition))
            SPositionArray.append(np.array(tempSPosition))
            PPositionArray.append(np.array(tempPPosition))
            # collection fluorescence
            # Donor channel
            ##Acceptor channel
            # fluoreAcceptor = fluorescence(dt, Qfluor = Qacceptor)
            # fluoreAcceptor.collectPhoton(molecularTrajectory.positionX, molecularTrajectory.positionY, molecularTrajectory.positionZ)


                # np.savetxt(f, fluoreDonor.singleTrace, fmt='%i')
            # with open('./acceptor_'+str(fileNum)+'.txt', 'a') as f:
            #    np.savetxt(f, fluoreAcceptor.trace, fmt='%i')
            #    #np.savetxt(f, fluoreAcceptor.singleTrace, fmt='%i')
            if t_cycle>t_flowin:
                allVolMolecules.append(volumeMolecule(molID,np.array([-border[0] / 2 + 0.5,
                                             np.random.uniform(-border[1] / 2, border[1] / 2),
                                            np.random.uniform(-border[2] / 2, border[2] / 2)]),states[np.random.binomial(1,pB)+1],None))
                molID+=1
                t_flowin+=dt_flowin
            t_cycle+=dt
        with open(path + '/moleculenum_' + str(fileNum) + '.txt', 'a') as f:
            np.savetxt(f, countsArray, fmt='%.3f')

    fig = plt.figure(figsize=(6, 6))
    Plot = axes3d.Axes3D(fig)#title="t="+str(0.001*totalTime*(repeatNum+fileNum*repeatCycle))+"ms")  # 创建一个三维的绘图工程
    Plot.view_init(25,-75)

    Plot.set_xlim3d([-10, 10])
    Plot.set_ylim3d([-2.5, 2.5])
    Plot.set_zlim3d([-0.5, 0.5])
    Plot.set_zlabel('Z')  # 坐标轴
    Plot.set_ylabel('Y')
    Plot.set_xlabel('X')

    PlotE=Plot.scatter(EnzymePositionArray[0][0], EnzymePositionArray[0][1], EnzymePositionArray[0][2], c='black')
    PlotS=Plot.scatter(SPositionArray[0][0], SPositionArray[0][1], SPositionArray[0][2], c='r')
    PlotP=Plot.scatter(PPositionArray[0][0], PPositionArray[0][1], PPositionArray[0][2], c='b')


    def update(n):
        global  EnzymePositionArray, SPositionArray, PPositionArray
        # n = frame counter

        # update the plot data
        PlotE.set_offsets(EnzymePositionArray[n])
        PlotS.set_offsets(SPositionArray[n])
        PlotP.set_offsets(PPositionArray[n])

    #  将数据点分成三部分画，在颜色上有区分度

    anim = animation.FuncAnimation(fig, update, frames=10,interval=10)
    anim.save(path+'/anitest.gif', writer='imagemagick')
    plt.savefig(path+'/screenshot.eps')
    plt.savefig(path + '/screenshot.png')
    plt.show()

with open(path + '/log.txt', 'a') as f:
    f.write(str(time.asctime()) + ' - trace: done.\n\n')
