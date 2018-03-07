import sys, os
import numpy as np
from truss import *

class NonLinearAlgorithm: # Father class for all the non-linear algorithms
    
    def __init__(self, _truss, _toll, _nItMax):
        self.truss = _truss                     # The truss
        self.toll = _toll                       # Tolerance used to assess convergence
        self.nItMax = _nItMax                   # Maximum number of iterations
        self.lamdaMax = 1.0                     # Maximum value for lamda (a priori, it can be different from 1)
        self.stopit = False                     # Only used by the GUI
        self.cleanWorkspace()                   # Run at the beginning of the simulation to clean the workspace: all the old result files will be deleted
        self.computeTangentMethod = 'analytic'  # Or 'numeric' -> Used in bar.buildKt()
        self.dhook = None                       # Only used by the GUI
        
    def display(self, step, lamda):
        if self.dhook:
            self.dhook.display(step, lamda)
            self.dhook.refresh()
        
    def archive(self, lamda, nIt): # Writes the files containing the results. Just the main information is archived, but you can enrich it if you want...
        f1 = open(('Lambda.ascii'),'a')
        f1.write(str(lamda)+'\n')
        f1.close()
        f2 = open(('Iterations.ascii'),'a')
        f2.write(str(nIt)+'\n')
        f2.close()
        for node in self.truss.nodes:
            f3 = open(('Node_'+str(node.nb)+'_POSITION.ascii'),'a')
            f3.write(str(node.x)+'   '+str(node.y)+'\n')
            f3.close()
            f4 = open(('Node_'+str(node.nb)+'_DISPLACEMENTS.ascii'),'a')
            f4.write(str(node.u)+'   '+str(node.v)+'\n')
            f4.close()
            f5 = open(('Node_'+str(node.nb)+'_F_EXT.ascii'),'a')
            f5.write(str(node.Fx)+'   '+str(node.Fy)+'\n')
            f5.close()
    
    def computeError(self): # Computes the error used to assess convergence
        error = 1.0
        # TODO: implement error computation here ...
        return error
    
    def cleanWorkspace(self):
        dir_name = os.getcwd()
        test = os.listdir(dir_name)
        
        for item in test:
            if item.endswith(".ascii"):
                os.remove(os.path.join(dir_name, item))
    
class IncrementalAlgorithm(NonLinearAlgorithm): # A simple explicit solver for non-linear problems
    
    def __init__(self, _truss, _dlamda): 
        
        NonLinearAlgorithm.__init__(self, _truss, 0., 0)
        self.dlamda = _dlamda
    
    def run(self):
        lamda = 0.
        lamda0 = 0.
        dlamda = self.dlamda
        step = 0
        while lamda < self.lamdaMax and not self.stopit:
            step+=1
            lamda = lamda0 + dlamda
            Delta_u = np.zeros((2*len(self.truss.nodes),1))
            
            print '\n--- Incremental method - Load level, lambda =', lamda, ' ---'
            
            self.truss.applyNodalLoads(lamda)
            
            Kt = self.truss.buildKt(self.computeTangentMethod)
            g = self.truss.getOOBF()
            du = np.linalg.solve(Kt, -g)
            
            Delta_u+=du
            self.truss.incrementPositions(Delta_u)
            
            self.truss.update()
            self.archive(lamda, 0)
            self.display(step, lamda)
            lamda0 = lamda

class NewtonRaphsonAlgorithm(NonLinearAlgorithm): # Newton-Raphson method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda): 
        
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax)
        self.dlamda = _dlamda
    
    def run(self):
        # TODO: implement Newton-Raphson method here ...
        return 1

class QuasiNewtonAlgorithm(NonLinearAlgorithm): # Quasi-Newton method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda): 
        
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax)
        self.dlamda = _dlamda
    
    def run(self):
        # TODO: implement Quasi-Newton method here ...
        return 1

class ArcLengthAlgorithm(NonLinearAlgorithm): # Arc-length method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id): 
        
        NonLinearAlgorithm.__init__(self, _truss, _toll, _nItMax)
        self.dlamda = _dlamda
        self.psi = _psi
        self.Id = _Id # Ideal number of iterations per step
        
    def run(self):
        # TODO: implement Arc-length method here ...
        return 1

class UpdatedNormalPlaneArcLengthAlgorithm(ArcLengthAlgorithm): # Updated normal plane arc-length method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id):
        ArcLengthAlgorithm.__init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id)
    
    def run(self):
        # TODO: implement Updated normal plane arc-length method here ...
        return 1

class NormalPlaneArcLengthAlgorithm(ArcLengthAlgorithm): # Normal plane arc-length method
    
    def __init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id):
        ArcLengthAlgorithm.__init__(self, _truss, _toll, _nItMax, _dlamda, _psi, _Id)
    
    def run(self):
        # TODO: implement Updated normal plane arc-length method here ...
        return 1