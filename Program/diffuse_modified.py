#simulation molecular trajectories and single-molecule fluorescence trace

from Program.trajectory import trajectory
from Program.fluorescence_wzq import fluorescence_wzq
from Program.reaction_3state import reaction_3state
import minpy.numpy as np
import numpy as np

import time
import sys

#paramaters initials
dt = 10 #time interval (us)

#simulation time (us)
#2**21: ~2 s
#2**20: ~1 s
#2**19: ~0.5 s
#2**18: ~0.25 s
#2**17: ~0.125 s
totalTime = 1000000

repeatCycle = 500 #repeat totalTime simulation
fileCycle = 3 #repeat simulation files

#simulatoin condition initialization
border = [3, 3, 15]
#moleculeNum =  np.int(sys.argv[3])*1000/12 #simulated molecule number
moleculeNum = 100

pA = np.float(sys.argv[11])
pB = np.float(sys.argv[12])

AmoleculeNum = int(moleculeNum*pA)
BmoleculeNum = int(moleculeNum*pB)
CmoleculeNum = moleculeNum - AmoleculeNum - BmoleculeNum

Adiffcoef = 2e-5*1125/np.int(sys.argv[7]) #component A diffusion coefficient (um^2/us), ~10ms
Bdiffcoef = 2e-5*1125/np.int(sys.argv[7]) #compent B diffusion coefficient (um^2/us)
Cdiffcoef = 2e-5*1125/np.int(sys.argv[7]) #compent B diffusion coefficient (um^2/us)
diffCoef = np.repeat([Adiffcoef, Bdiffcoef, Cdiffcoef], [AmoleculeNum, BmoleculeNum, CmoleculeNum])

#fluor quantum yield
Qfluor=1
QA = np.float(sys.argv[8])
QB = np.float(sys.argv[9])
QC = np.float(sys.argv[10])
QacceptorA = 0
QdonorB = 0
QacceptorB = 0

#reaction detail
#Input: k01 k10 k12 k21 k20 k02
k_Matrix=np.zeros((3,3))
k_Matrix[0][1]=np.int(sys.argv[1])/1e6 #us-1
k_Matrix[1][0]=np.int(sys.argv[2])/1e6 #us-1
k_Matrix[1][2]=np.int(sys.argv[3])/1e6 #us-1
k_Matrix[2][1]=np.int(sys.argv[4])/1e6 #us-1
k_Matrix[2][0]=np.int(sys.argv[5])/1e6 #us-1
k_Matrix[0][2]=np.int(sys.argv[6])/1e6 #us-1
k_Matrix[0][0]=-k_Matrix[0][1]-k_Matrix[0][2]
k_Matrix[1][1]=-k_Matrix[1][0]-k_Matrix[1][2]
k_Matrix[2][2]=-k_Matrix[2][0]-k_Matrix[2][1]

initInfo =  'trajectory simulation time interval: '+str(dt)+' us\n'+ \
            'total simulation time: '+str(totalTime*1e-6)+'x'+str(repeatCycle)+' s\n'+ \
            'molecule number: '+str(moleculeNum)+'\n'+ \
            'border volume: '+'x'.join(map(str, border))+' um^3\n'+ \
            'apparent concentration: '+str(int(moleculeNum/6.023/border[0]/border[1]/border[2]*1e4))+' pM\n'+ \
            '-------------------------------------\n'+ \
            'A ratio: '+str(pA)+'\n'+ \
            'B ratio: '+str(pB)+'\n'+ \
            'A molecule number: '+str(AmoleculeNum)+'\n'+ \
            'B molecule number: '+str(BmoleculeNum)+'\n'+ \
            'C molecule number: '+str(CmoleculeNum)+'\n'+ \
            'A diffusion coefficient: '+str(Adiffcoef*1e-6) +' m^2/s\n'+ \
            'B diffusion coefficient: '+str(Bdiffcoef*1e-6) +' m^2/s\n'+ \
            'C diffusion coefficient: '+str(Cdiffcoef*1e-6) +' m^2/s\n'+ \
            'A effective brightness: '+str(Qfluor*4.172)+'\n'+ \
            '-------------------------------------\n\n'+ \
            'Reaction rate constant : '+(' ').join(map(lambda x: (' ').join(map(str,x)),k_Matrix))+'\n'+ \
            'Intensity A:'+str(QA)+'\n'+ \
            'Intensity B:'+str(QB)+'\n'+ \
            '-------------------------------------\n\n'+ \
            '1-typical, 2-without poisson random, 3-typical without int, 4-without poisson random & int' + \
            'without poisson noise & crosstalk.\n\n'

#log information file;

path = '../180421serie/D'+sys.argv[7]+'_QA'+str(QA)+'_QB'+str(QB)+'_QC'+str(QC)+'_pA'+str(pA)+'_pB'+str(pB)

with open(path + '/log.txt', 'w') as f:
    f.write(initInfo+str(time.asctime())+' - trace: simulating...\n\n')

for fileNum in range(fileCycle):
    #initial stochastic coordinates generation
    initPosition = [np.random.uniform(-border[0]/2, border[0]/2, moleculeNum),
                    np.random.uniform(-border[1]/2, border[1]/2, moleculeNum),
                    np.random.uniform(-border[2]/2, border[2]/2, moleculeNum)]
    initState=np.zeros(moleculeNum)
    for i in range (moleculeNum): initState[i]= QA if i<AmoleculeNum else(QB if (i>=AmoleculeNum and i<AmoleculeNum+BmoleculeNum) else QC)
       
    #donor channel trace file;
    #accptor channel trace file;
    with open(path + '/donor_'+str(fileNum)+'.txt', 'a') as f:
        f.write('')
            
    with open(path + '/donor_'+str(fileNum)+'nr.txt', 'a') as f:
        f.write('')
   
    for repeatNum in range(repeatCycle):
        print('Simulating: Cycle'+str(repeatNum)+"\n")
        #trajectory simulation
        molecularTrajectory = trajectory(dt, totalTime, moleculeNum, initPosition, border = border)
        molecularTrajectory.diffusion(diffCoef)
        initPosition = [molecularTrajectory.positionX[-1, :],
                        molecularTrajectory.positionY[-1, :],
                        molecularTrajectory.positionZ[-1, :]]

        #reaction simulation
        if k_Matrix[0][1] > 0:
            reactionTrajectory=reaction_3state(dt, totalTime, moleculeNum, initState)
            reactionTrajectory.react(k_Matrix,QA,QB,QC)
            initState=reactionTrajectory.state[-1, :]
            fluoreDonor = fluorescence_wzq(dt, Qfluor)
            fluoreDonor.collectPhoton(molecularTrajectory.positionX, molecularTrajectory.positionY, molecularTrajectory.positionZ,reactionTrajectory.state)
        else:
            state=np.ones([totalTime/dt,moleculeNum])
            fluoreDonor = fluorescence_wzq(dt)
            fluoreDonor.collectPhoton(molecularTrajectory.positionX, molecularTrajectory.positionY, molecularTrajectory.positionZ,state)
                
        #collection fluorescence
        #Donor channel
        ##Acceptor channel
        #fluoreAcceptor = fluorescence(dt, Qfluor = Qacceptor)
        #fluoreAcceptor.collectPhoton(molecularTrajectory.positionX, molecularTrajectory.positionY, molecularTrajectory.positionZ)

        with open(path + '/donor_'+str(fileNum)+'.txt', 'a') as f:
            np.savetxt(f, fluoreDonor.trace[:-10000], fmt='%.3f')
            
        with open(path + '/donor_'+str(fileNum)+'nr.txt', 'a') as f:
            np.savetxt(f, fluoreDonor.trace_nr[:-10000], fmt='%.3f')

        with open(path + '/moleculenum_' + str(fileNum) + '.txt', 'a') as f:
            np.savetxt(f, reactionTrajectory.moleculenum[:-10000], fmt='%.3f')
            #np.savetxt(f, fluoreDonor.singleTrace, fmt='%i')
        #with open('./acceptor_'+str(fileNum)+'.txt', 'a') as f:
        #    np.savetxt(f, fluoreAcceptor.trace, fmt='%i')
        #    #np.savetxt(f, fluoreAcceptor.singleTrace, fmt='%i')
    
with open(path + '/log.txt', 'a') as f:
    f.write(str(time.asctime())+' - trace: done.\n\n')
