from node import *
from bar import *
import numpy as np

class Truss():

    def __init__(self):
        self.nodes = []
        self.bars = []
        self.fixedNodes = []
    
    def addNode(self, node): # Adds a node to the truss
        self.nodes.append(node)
    
    def addBar(self, bar): # Adds a bar to the truss
        self.bars.append(bar)
    
    def setNodeRows(self): # Set the rows corresponding to the nodes degrees of freedom in the global system of equations -> Used to imposed boundary conditions
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            node.rowX = 2*i
            node.rowY = 2*i+1
        
    def buildKt(self, computeTangentMethod): # Assembles the global tangent stiffness matrix
        Kt = np.zeros((2*len(self.nodes), 2*len(self.nodes)))
        for bar in self.bars:
            Kte = bar.buildKt(computeTangentMethod)
            for i in range(len(bar.nodes)):
                nodi = bar.nodes[i]
                for j in range(len(bar.nodes)):
                    nodj = bar.nodes[j]
                    
                    Kt[nodi.rowX, nodj.rowX] += Kte[2*i, 2*j]
                    Kt[nodi.rowX, nodj.rowY] += Kte[2*i, 2*j+1]
                    Kt[nodi.rowY, nodj.rowX] += Kte[2*i+1, 2*j]
                    Kt[nodi.rowY, nodj.rowY] += Kte[2*i+1, 2*j+1]
        
        self.applyBCs(Kt)
        return Kt
    
    def getFint(self): # Assembles the global internal forces vector
        Fint = np.zeros((2*len(self.nodes), 1))
        for bar in self.bars:
            Fint_e = bar.getFint()
            for i in range(len(bar.nodes)):
                node = bar.nodes[i]
                Fint[node.rowX]+=Fint_e[2*i]
                Fint[node.rowY]+=Fint_e[2*i+1]
        return Fint
        
    def getFext(self): # Assembles the global external forces vector
        Fext = np.zeros((2*len(self.nodes), 1))
        for node in self.nodes:
            Fext[node.rowX] = node.Fx
            Fext[node.rowY] = node.Fy
        return Fext
        
    def get_qef(self): # Assembles the global loads vector
        qef = np.zeros((2*len(self.nodes), 1))
        for node in self.nodes:
            qef[node.rowX] = node.qef_x
            qef[node.rowY] = node.qef_y
        return qef
    
    def getOOBF(self): # Computes the Out-Of-Balance-Forces vector
        g = np.zeros((2*len(self.nodes), 1))
        
        Fint = self.getFint()
        Fext = self.getFext()
        g = Fint - Fext
        
        for node in self.fixedNodes:
            if node.isFixedAlongX:
                g[node.rowX] = 0.
            if node.isFixedAlongY:
                g[node.rowY] = 0.
        return g
    
    def fix(self, node, dof): # Imposes boundary conditions (only fixations are considered for the moment)
        self.fixedNodes.append(node)
        if (dof == 'x'):
            node.isFixedAlongX = True
        elif (dof == 'y'):
            node.isFixedAlongY = True
        else:
            raise Exception('Unknown dof!')
        
    def applyBCs(self, Kt): # Applies boundary conditions on the global tangent stiffness matrix
        for node in self.fixedNodes:
            if node.isFixedAlongX:
                Kt[node.rowX, :] = 0.
                Kt[:, node.rowX] = 0.
                Kt[node.rowX, node.rowX] = 1.0
            if node.isFixedAlongY:
                Kt[node.rowY, :] = 0.
                Kt[:, node.rowY] = 0.
                Kt[node.rowY, node.rowY] = 1.0
    
    def applyNodalLoads(self, lamda): # Applies the current loads to the nodes
        for node in self.nodes:
            node.Fx = node.qef_x * lamda
            node.Fy = node.qef_y * lamda
    
    def applyNodalTimeLoads(self, t): # Applies the current time-dependent loads to the nodes --> NB: To use only with NewtonRaphsonAlgorithmInTime !
        for node in self.nodes:
            node.Fx = node.get_qef_x(t)
            node.Fy = node.get_qef_y(t)
    
    def incrementPositions(self, du): # Increments nodes displacements of a quantity 'du' (which is a vector containing all the displacement variations) and updates nodes positions according to the new displacements
        for node in self.nodes:
            node.u = node.uOld + du[node.rowX,0]
            node.v = node.vOld + du[node.rowY,0]
            node.x = node.x0 + node.u
            node.y = node.y0 + node.v
    
    def resetPositions(self): # Reset nodes displacements and positions to the ones at the beginning of the step -> Useful when you have to restart a step!
        for node in self.nodes:
            node.u = node.uOld
            node.v = node.vOld
            node.x = node.xOld
            node.y = node.yOld
    
    def resetNodalLoads(self): # Resets current nodes loads to zero
        for node in self.nodes:
            node.Fx = 0.
            node.Fy = 0.
    
    def update(self):
        for node in self.nodes: # Updates nodes quantities -> To be called at the end of a step (i.e. when convergence is reached)
            node.uOld = node.u
            node.vOld = node.v
            node.xOld = node.x
            node.yOld = node.y