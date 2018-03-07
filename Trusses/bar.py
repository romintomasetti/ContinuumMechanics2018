from node import *
import numpy as np
import math

class Bar:
    
    def __init__(self, _nb, nods, _E, _A):
        self.nb = _nb                                                   # Bar number
        self.nodes = nods                                               # Bar nodes (it is an array)
        if not len(self.nodes) == 2:
            raise Exception('bar nodes number is different from two!')
        self.l0 = self.getLength()                                      # Bar initial length
        self.E = _E                                                     # Bar Young modulus
        self.A = _A                                                     # Bar section area (always supposed constant in this code)
    
    def getLength(self):
        return math.sqrt((self.nodes[1].x - self.nodes[0].x)**2 + (self.nodes[1].y - self.nodes[0].y)**2)
    
    def getE_GL(self):
        E_GL = 0. 
        # TODO: implement the computation of Gree-Lagrange strains here ...
        
        raise Exception('Implement the computation of Green-Lagrange strains!') # <-- To be deleted once the computation of Green-Lagrange strains implemented!
        
        return E_GL
        
    def getPK2Stress(self):
        PK2 = 0. 
        # TODO: implement the computation of PK2 stresses here ...
        
        raise Exception('Implement the computation of PK2 stresses!') # <-- To be deleted once the computation of PK2 stresses implemented!
        
        return PK2
    
    def getFint(self): # Computes the bar internal forces vector
        
        Fint = np.zeros((0, 0)) # --> Replace the zeros with the right vector size!
        
        # TODO: implement the computation of the internal forces here ...
        
        raise Exception('Implement the computation of the internal forces!') # <-- To be deleted once the computation of internal forces implemented!
        
        return Fint
    
    def buildKt(self, computeTangentMethod): # Computes the bar tangent stiffness matrix
        
        Kt = np.zeros((0, 0)) # --> Replace the zeros with the right matrix size!
        
        # TODO: implement the computation of the tangent stiffness matrix here ...
        
        raise Exception('Implement the computation of the tangent stiffness matrix!') # <-- To be deleted once the computation of tangent stiffness matrix implemented!
        
        return Kt