import numpy as np

class state:
    def __init__(self,ID,brightness,diffconstant):
        self.ID=ID
        self.brightness=brightness
        self.Topo=np.array([])
        self.diffconstant=diffconstant
    def chemlink(self,otherstate):
        self.Topo = np.append(self.Topo,otherstate)

class molecule:
    def __init__(self,ID,initposition,initstate,initbinding):
        self.ID=ID
        self.position=np.zeros(3)
        self.position=initposition
        self.state=initstate
        self.bindingmolecule=initbinding
        self.uniTransitTime=None
        self.toState=None
        self.time=0
    def getID(self):
        return self.ID
    def getposition(self):
        return self.position
    def getstate(self):
        return self.state

class surfMolecule(molecule):
    def __init__(self,ID,initposition,initstate,initbinding,rigradius,cbradius,ubradius,k_Matrix):
        molecule.__init__(self,ID,initposition,initstate,initbinding)
        self.rigradius=rigradius
        self.cbradius=cbradius
        self.rRelax=ubradius
        self.k_Matrix=k_Matrix
    def unireactTest(self):
        if self.uniTransitTime==None and self.state.ID!=0:
            waitingTime = np.random.exponential(1 / (-self.k_Matrix[self.state.ID][self.state.ID]))
            if waitingTime > 0:
                self.uniTransitTime= self.time + waitingTime
            jump_judge = self.k_Matrix[self.state.ID][self.state.Topo[0].ID] / (-self.k_Matrix[self.state.ID][self.state.ID])
            if np.random.uniform(0, 1) < jump_judge:
                self.toState=self.state.Topo[0]
            else:
                self.toState=self.state.Topo[1]
        return self.uniTransitTime

    def bind(self, TargetMolecule:molecule):
        self.bindingmolecule=TargetMolecule
        TargetMolecule.bind(self)
        self.state=TargetMolecule.state

    def bireactTest(self,Targetmolecule:molecule,dt):
        if self.state.ID==0 and Targetmolecule.isactive:
            r=np.linalg.norm(Targetmolecule.position-self.position)
            if r>self.rigradius: pr=dt*self.k_Matrix[self.state.ID][Targetmolecule.state.ID]/(2*np.pi*self.cbradius**3)*np.exp(-(r-self.rigradius)/self.cbradius)
            elif r<=self.rigradius: pr=dt*self.k_Matrix[self.state.ID][Targetmolecule.state.ID]/(2*np.pi*self.cbradius**3)
            if np.random.uniform()<=pr:
                self.bind(Targetmolecule)

    def unbind(self):
        self.bindingmolecule.unbind()
        self.bindingmolecule=None

    def changeState(self):
        if self.toState!=None:
            if self.toState.ID==0:
                self.unbind()
                self.state=self.toState
            elif self.toState.ID==1 or self.toState.ID==2:
                self.state=self.toState
                self.bindingmolecule.state=self.toState
            self.toState=None
            self.uniTransitTime=None

    def update(self,dt):
        self.time+=dt
        if self.uniTransitTime:
            if self.time>=self.uniTransitTime:
                self.changeState()
        self.unireactTest()
        return True

class volumeMolecule(molecule):
    def __init__(self,ID,initposition,initstate,initbinding):
        molecule.__init__(self,ID,initposition,initstate,initbinding)
        self.isactive=True
        self.inbox=True

    def resetState(self,state):
        self.state=state

    def bind(self, targetMolecule:molecule):
        self.bindingmolecule=targetMolecule
        self.isactive=False
        self.position[0] = targetMolecule.position[0]
        self.position[1] = targetMolecule.position[1]
        self.position[2] = targetMolecule.position[2]

    def unbind(self):
        theta=np.random.uniform(0,2*np.pi)
        phi=np.random.uniform(0,np.pi)
        self.position[0] = self.bindingmolecule.rRelax*np.cos(theta)*np.sin(phi)
        self.position[1] = self.bindingmolecule.rRelax*np.sin(theta)*np.sin(phi)
        self.position[2] = self.bindingmolecule.rRelax*np.cos(phi)
        self.bindingmolecule=None
        self.isactive=True


    def diffuse(self,dt,border,drift):
        if self.isactive:
            walkdist=np.sqrt(2*self.state.diffconstant*dt)
            self.position=np.random.normal(self.position+drift*dt,walkdist)
        if np.abs(self.position[1]) > border[1]/2:
            self.position[1]=self.position[1]/np.abs(self.position)[1]*border[1]-self.position[1]
        if np.abs(self.position[2]) > border[2]/2:
            self.position[2]=self.position[2]/np.abs(self.position)[2]*border[2]-self.position[2]
        if self.position[0] < -border[0]/2:
            self.position[0] = -border[0]-self.position[0]
        elif np.abs(self.position[0]) > border[0]/2:
            self.position[0] = -border[0] + self.position[0]
            self.inbox=False

    def isinbox(self):
        return self.inbox

    def update(self,dt,extborder,extdrift):
        self.time+=dt
        self.diffuse(dt,border=extborder,drift=extdrift)
