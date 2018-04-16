#simulation molecular trajectories and single-molecule fluorescence trace

from trajectory import trajectory
from fluorescence_wzq import fluorescence_wzq
from reaction_wzq import reaction_wzq
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

repeatCycle = 2000 #repeat totalTime simulation
fileCycle = 3 #repeat simulation files

#simulatoin condition initialization
border = [3, 3, 15]
#moleculeNum =  np.int(sys.argv[3])*1000/12 #simulated molecule number
moleculeNum = 30

pA = 1
AmoleculeNum = int(moleculeNum*pA)
BmoleculeNum = moleculeNum - AmoleculeNum

Adiffcoef = 2e-5*1125/np.int(sys.argv[4]) #compent A diffusion coefficient (um^2/us), ~10ms
Bdiffcoef = 6e-5 #compent B diffusion coefficient (um^2/us)
diffCoef = np.repeat([Adiffcoef, Bdiffcoef], [AmoleculeNum, BmoleculeNum])

#fluor quantum yield
QA = np.float(sys.argv[5])
Q = np.float(sys.argv[3]) # = light intensity of dark component/light intensity of bright component
QacceptorA = 0
QdonorB = 0
QacceptorB = 0

#reaction detail
kplus=np.int(sys.argv[1])/1e6 #us-1
kminus=np.int(sys.argv[2])/1e6 #us-1

initInfo =  'trajectory simulation time interval: '+str(dt)+' us\n'+ \
            'total simulation time: '+str(totalTime*1e-6)+'x'+str(repeatCycle)+' s\n'+ \
            'molecule number: '+str(moleculeNum)+'\n'+ \
            'border volume: '+'x'.join(map(str, border))+' um^3\n'+ \
            'apparent concentration: '+str(int(moleculeNum/6.023/border[0]/border[1]/border[2]*1e4))+' pM\n'+ \
            '-------------------------------------\n'+ \
            'A ratio: '+str(pA)+'\n'+ \
            'A molecule number: '+str(AmoleculeNum)+'\n'+ \
            'B molecule number: '+str(BmoleculeNum)+'\n'+ \
            'A diffusion coefficient: '+str(Adiffcoef*1e-6) +' m^2/s\n'+ \
            'B diffusion coefficient: '+str(Bdiffcoef*1e-6) +' m^2/s\n'+ \
            'A effective brightness: '+str(QA*4.172)+'\n'+ \
            '-------------------------------------\n\n'+ \
            'forward reaction rate constant (bright to dark): '+str(kplus)+'\n'+ \
            'reverse reaction rate constant (dark to bright) '+str(kminus)+'\n'+ \
            'relative intensity '+str(Q)+'\n'+ \
            '-------------------------------------\n\n'+ \
            '1-typical, 2-without poisson random, 3-typical without int, 4-without poisson random & int' + \
            'without poisson noise & crosstalk.\n\n'

#log information file;

path = './160518serie/D'+sys.argv[4]+'_kplus'+str(np.int(kplus*1e6))+'_kminus'+str(np.int(kminus*1e6))+'_Q'+str(Q)+'_QA'+str(QA)

with open(path + '/log.txt', 'w') as f:
    f.write(initInfo+str(time.asctime())+' - trace: simulating...\n\n')

for fileNum in range(fileCycle):
    #initial stochastic coordinates generation
    initPosition = [np.random.uniform(-border[0]/2, border[0]/2, moleculeNum),
                    np.random.uniform(-border[1]/2, border[1]/2, moleculeNum),
                    np.random.uniform(-border[2]/2, border[2]/2, moleculeNum)]
       
    #donor channel trace file;
    #accptor channel trace file;
    with open(path + '/donor_'+str(fileNum)+'.txt', 'a') as f:
        f.write('')
            
    with open(path + '/donor_'+str(fileNum)+'nr.txt', 'a') as f:
        f.write('')
   
    for repeatNum in range(repeatCycle):
        #trajectory simulation
        molecularTrajectory = trajectory(dt, totalTime, moleculeNum, initPosition, border = border)
        molecularTrajectory.diffusion(diffCoef)
        initPosition = [molecularTrajectory.positionX[-1, :],
                        molecularTrajectory.positionY[-1, :],
                        molecularTrajectory.positionZ[-1, :]]

        #reaction simulation
        if kplus > 0:
            reactionTrajectory=reaction_wzq(dt, totalTime, moleculeNum)
            reactionTrajectory.react(kplus,kminus,Q)
            fluoreDonor = fluorescence_wzq(dt, Qfluor=QA)
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
            np.savetxt(f, fluoreDonor.trace, fmt='%.3f')
            
        with open(path + '/donor_'+str(fileNum)+'nr.txt', 'a') as f:
            np.savetxt(f, fluoreDonor.trace_nr, fmt='%.3f')
            #np.savetxt(f, fluoreDonor.singleTrace, fmt='%i')
        #with open('./acceptor_'+str(fileNum)+'.txt', 'a') as f:
        #    np.savetxt(f, fluoreAcceptor.trace, fmt='%i')
        #    #np.savetxt(f, fluoreAcceptor.singleTrace, fmt='%i')
    
with open(path + '/log.txt', 'a') as f:
    f.write(str(time.asctime())+' - trace: done.\n\n')
