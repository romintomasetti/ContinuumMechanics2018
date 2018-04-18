from node import *
import numpy as np
import math

class Bar:
    
    def __init__(self, _nb, nods, _E, _A):
        self.nb = _nb                                                   # Bar number
        self.nodes = nods                                               # Bar nodes (it is an array)
        if not len(self.nodes) == 2:
            raise Exception('bar nodes number is different from two!')
        self.l0        = self.getLength()                               # Bar initial length
        self.E_modulus = _E                                             # Bar Young modulus
        self.Area      = _A                                             # Bar section area (always supposed constant in this code)
        self.numberOfKtEvals = 0   
    
    def getLength(self):
        return math.sqrt((self.nodes[1].x - self.nodes[0].x)**2 + (self.nodes[1].y - self.nodes[0].y)**2)
    
    def getE_GL(self):
        E_GL = ( self.getLength()**2 - self.l0**2 ) / ( 2 * self.l0**2 )
        return E_GL
        
    def getPK2Stress(self): # Computes PK2 stress
        E_GL = self.getE_GL()
        PK2  = self.E_modulus * E_GL
        #print " -------------- Pk2 = ", self.E_modulus*E_GL, "--------------"
        return PK2
    
    def getFint(self): # Computes the bar internal forces vector
        
        Fint       = np.zeros(4) 
        PK2_stress = self.getPK2Stress()
        C          = self.get_C()
        constant   = PK2_stress * self.Area

        Fint[0] = -C[0]
        Fint[1] = -C[1]
        Fint[2] =  C[0]
        Fint[3] =  C[1]
        
        Fint = constant * Fint
        
        return Fint
    
    def buildKt(self, computeTangentMethod): # Computes the bar tangent stiffness matrix

        if computeTangentMethod == 'numeric':
            self.numberOfKtEvals += 1
            
            delta = 1e-5
            
            Kt = np.zeros((4,4))
            
            for j in range(0,4):
                if j >= 2:
                    n = 1;
                else:
                    n = 0;
                if j%2 == 0:
                    self.nodes[n].x += delta
                else:
                    self.nodes[n].y += delta
                P_plus = self.getFint()
                if j%2 == 0:
                    self.nodes[n].x -= 2*delta
                else:
                    self.nodes[n].y -= 2*delta
                P_minus = self.getFint()
                
                print P_plus,P_minus
                
                if j%2 == 0:
                    self.nodes[n].x += delta
                else:
                    self.nodes[n].y += delta
                for i in range(0,4):
                    Kt[i][j] = (P_plus[i]-P_minus[i])/(2*delta)
                    
            print Kt
            
        else:
            self.numberOfKtEvals += 1
            # The tangent stiffness matrix of a bar is 4 by 4:
            Kt = np.zeros((4, 4))
            
            # The tangent stiffness matrix is a combination of two matrices.
            # The material stiffness matrix:
            K_mat = np.zeros((4,4))
            # The geometric stiffness matrix:
            K_geo = np.zeros((4,4))
            
            # Building of K_mat (see https://www.colorado.edu/engineering/cas/courses.d/NFEM.d/NFEM.Ch09.d/NFEM.Ch09.pdf)
            constant_mat = self.Area * self.E_modulus / self.l0
            C = self.get_C()
            
            K_mat[0][0] =  C[0]**2
            
            K_mat[0][1] =  C[0] * C[1]
            K_mat[1][0] =  K_mat[0][1]
            
            K_mat[0][2] = -C[0]**2
            K_mat[2][0] =  K_mat[0][2]
            
            K_mat[0][3] = -C[0] * C[1]
            K_mat[3][0] =  K_mat[0][3]
            
            K_mat[1][1] =  C[1]**2
            
            K_mat[1][2] =  K_mat[0][3]
            K_mat[2][1] =  K_mat[1][2]
            
            K_mat[1][3] = -C[1]**2
            K_mat[3][1] = K_mat[1][3]
            
            K_mat[2][2] = K_mat[0][0]
            
            K_mat[2][3] = K_mat[0][1]
            K_mat[3][2] = K_mat[2][3]
            
            K_mat[3][3] = K_mat[1][1]
            
            K_mat = constant_mat * K_mat
            
            # Building of K_geo:
            # Get PK2 stress:
            PK2_stress = self.getPK2Stress()
            constant_geo = self.Area * PK2_stress / self.l0
            
            K_geo[0][0] = 1
            K_geo[1][1] = 1
            K_geo[2][2] = 1
            K_geo[3][3] = 1
            K_geo[0][2] = -1
            K_geo[1][3] = -1
            K_geo[2][0] = -1
            K_geo[3][1] = -1
            
            K_geo = constant_geo * K_geo
            
            Kt = K_geo + K_mat        
        
        return Kt
        
    def get_C(self):
        C = np.zeros(2)
        C[0] = (self.nodes[1].x-self.nodes[0].x)/self.l0
        C[1] = (self.nodes[1].y-self.nodes[0].y)/self.l0
        return C
        
   
            
        