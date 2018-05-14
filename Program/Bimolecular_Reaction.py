import numpy as np

class state:
    def __init__(self,ID,brightness,diffconstant):
        self.ID=ID
        self.brightness=brightness
        self.Topo=np.array([])
        self.diffconstant=diffconstant
    def chemlink(self,otherstate):
        np.append(self.Topo,otherstate)

class molecule:
    def __init__(self,initposition,initstate,initbinding):
        self.position=np.zeros(3)
        self.position=initposition
        self.state=initstate
        self.bindingmolecule=initbinding
        self.uniTransitTime=None
        self.toState=None
        self.time=0

class surfMolecule(molecule):
    def __init__(self,initposition,initstate,initbinding,cbradius,k_Matrix):
        molecule.__init__(self,initposition,initstate,initbinding)
        self.cbradius=cbradius
        self.k_Matrix=k_Matrix
    def unireactTest(self):
        if self.uniTransitTime==None and self.state.ID!=0:
            waitingTime = np.random.exponential(1 / (-self.k_Matrix[self.state.ID][self.state.ID]))
            if waitingTime > 0:
                self.uniTransitTime= self.time + waitingTime
            jump_judge = self.k_Matrix[self.state.ID][self.state.Topo[0]] / (-self.k_Matrix[self.state.ID][self.state.ID])
            if np.random.uniform(0, 1) < jump_judge:
                self.toState=self.state.Topo[0]
            else:
                self.toState=self.state.Topo[1]
        return self.uniTransitTime

    def bireactTest(self,Targetmolecule:molecule,dt):
        if self.state.ID==0 and Targetmolecule.isactive:
            r=np.linalg.norm(Targetmolecule.position-self.position)
            pr=dt*self.k_Matrix[self.state.ID][Targetmolecule.state.ID]/(2*np.pi*self.cbradius**3)*np.exp(-r/self.cbradius)
            if np.random.uniform()<=pr:
                self.bind(self,Targetmolecule)

    def bind(self, TargetMolecule:molecule):
        self.bindingmolecule=TargetMolecule
        TargetMolecule.bind(self)
        self.state=TargetMolecule.state

    def unbind(self):
        self.bindingmolecule.bindingmolecule=None
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
        if self.time>=self.uniTransitTime:
            self.changeState()
        self.unireactTest()
        return True

class volumeMolecule(molecule):
    def __init__(self,initposition,initstate,initbinding):
        molecule.__init__(self,initposition,initstate,initbinding)
        self.positionmemory=None
        self.isactive=True
        self.velocity=np.zeros(3)

    def resetState(self,state):
        self.state=state

    def bind(self, targetMolecule:molecule):
        self.bindingmolecule=targetMolecule
        self.isactive=False
        self.positionmemory=self.position
        self.position=targetMolecule.position

    def unbind(self):
        self.bindingmolecule=None
        self.isactive=True
        self.position=self.positionmemory
        self.positionmemory=None

    def diffuse(self,dt,border):
        if self.isactive:
            walkdist=np.sqrt(2*self.state.diffconstant*dt)
            self.position=np.random.normal(self.position+self.velocity*dt,walkdist)
        ind=np.abs(self.position) > border/2
        if np.sum(ind):
            self.position[ind]=self.position[ind]/np.abs(self.position)[ind]*border[ind]-self.position[ind]
    def update(self,dt):
        self.time+=dt
        self.diffuse(dt)
