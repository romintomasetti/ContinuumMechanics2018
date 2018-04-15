class Node:
    
    def __init__(self, _nb, _x, _y):
        self.nb = _nb               # Node number -> NB: Each node has to have a different number!
        self.x = _x                 # Node current position along x
        self.y = _y                 # Node current position along y
        self.xOld = self.x          # Node position along x at the beginning of the step
        self.yOld = self.y          # Node position along y at the beginning of the step
        self.x0 = self.x            # Node position along x at the beginning of the simulation ("t=0")
        self.y0 = self.y            # Node position along y at the beginning of the simulation ("t=0")
        self.u = 0.                 # Node current displacement along x
        self.v = 0.                 # Node current displacement along y
        self.uOld = self.u          # Node displacement along x at the beginning of the step
        self.vOld = self.v          # Node displacement along y at the beginning of the step
        self.rowX = -1              # Row corresponding to the node x degree of freedom in the global system of equations -> Used to imposed boundary conditions, set using truss.setNodeRows()
        self.rowY = -1              # Row corresponding to the node y degree of freedom in the global system of equations -> Used to imposed boundary conditions, set using truss.setNodeRows()
        self.Fx = 0.                # Load currently applied on the node along x (= qef_x * lamda)
        self.Fy = 0.                # Load currently applied on the node along y (= qef_x * lamda)
        self.qef_x = 0.             # Total (initial) load applied on the node along x
        self.qef_y = 0.             # Total (initial) load applied on the node along y
        self.qef_x_t = None         # Time-dependent load applied on the node along x (it's a python function)
        self.qef_y_t = None         # Time-dependent load applied on the node along y (it's a python function)
        self.isFixedAlongX = False  # Flag set to true if the node is fixed along x
        self.isFixedAlongY = False  # Flag set to true if the node is fixed along y
    
    def applyLoad(self, dof, val):  # Applies a load of intensity 'val' along the direction identified by 'dof'
        if dof == 'x':
            self.qef_x = val
        if dof == 'y':
            self.qef_y = val
    
    def applyTimeLoad(self, dof, _f): # Initialize the time-dependent load (which is simply a python function) along the direction identified by 'dof' --> NB: To use only with NewtonRaphsonAlgorithmInTime !
        if dof == 'x':
            self.qef_x_t = _f
        if dof == 'y':
            self.qef_y_t = _f
    
    def get_qef_x(self, t): # Evaluates the time-dependent load along x at time 't' --> NB: To use only with NewtonRaphsonAlgorithmInTime !
        if self.qef_x_t:
            return self.qef_x_t(t)
        else:
            return 0.
    
    def get_qef_y(self, t): # Evaluates the time-dependent load along y at time 't' --> NB: To use only with NewtonRaphsonAlgorithmInTime !
        if self.qef_y_t:
            return self.qef_y_t(t)
        else:
            return 0.