import numpy as np

class state:
    def __init__(self,ID,brightness,chemtopo,kineticmatrix):
        self.ID=ID
        self.brightness=brightness


class molecule:
    def __init__(self,initposition,initstate,initbinding,diffconstant):
        self.position=np.zeros(3)
        self.position=initposition
        self.state=initstate
        self.bindingmolecule=initbinding
        self.transitTime=None
        self.toState=None
        self.time=0
        self.diffconstant=diffconstant
class surfmolecule(molecule):
    def __init__(self,initposition,initstate,initbinding):
        molecule.__init__(self,initposition,initstate,initbinding)
    def unireact(self):

    def bireact(self,bindingmolecule):

    def bind(self,bindingmolecule):
        self.bindingmolecule=bindingmolecule
        bindingmolecule.bindingmolecule=self
        self.state=bindingmolecule.state
    def unbind(self):
        self.bindingmolecule.bindingmolecule=None
        self.bindingmolecule.position=np.random.normal()
        self.bindingmolecule=None
    def changeState(self):
        if self.toState.ID==0:
            self.unbinding()
            self.state=self.toState
