import sys
sys.path.insert(0, '../../')
from node import *
from bar import *
from truss import *
from nonLinAlgo import *

import matplotlib.pyplot as plt
import math

from optparse import OptionParser
import vtk

def main(nogui):
    
    #Geometry
    a = 2.0
    b = 1.0
    #sys.exit()
    l0 = math.sqrt(a**2 + b**2)
    tol = 10**(-6)
    nItMax = 100
    
    #Nodes
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, a, b)
    #node3 = Node(3, 2*a, 0.0)
    
    #Material properties
    E = 70e9
    A = 0.01
    
    #Critical load
    qcr = (math.sqrt(3)*E*A*(b**3))/(9.0*(l0**3)) #It's the critical load divided by 2 since we consider just one bar!
    
    #Bars
    bar1 = Bar(1, [node1, node2], E, A)
    #bar2 = Bar(2, [node2, node3], E, A)
    
    #Truss
    truss = Truss()
    truss.addNode(node1)
    truss.addNode(node2)
    #truss.addNode(node3)
    truss.addBar(bar1)
    #truss.addBar(bar2)
    truss.setNodeRows() # ATTENTION: this line is mandatory at the end of the definition of the truss!
    
    #BCs
    truss.fix(node1,'x')
    truss.fix(node1,'y')
    truss.fix(node2,'x')
    #truss.fix(node3,'x')
    #truss.fix(node3,'y')
    
    #Loads
    node2.applyLoad('y', -1.5*qcr)
    
    #Non-linear algorithm
    dlamda0 = 0.01
    algo = IncrementalAlgorithm(truss, dlamda0)
    
    #GUI --> if you want to run a test without GUI: in a command window type (without the quotation marks) 'python test.py --nogui'
    if nogui:
        if vtk.VTK_MAJOR_VERSION <= 5:
            import trussViewerVtk5PyQt4 as v
        elif vtk.VTK_MAJOR_VERSION == 6:
            import trussViewerVtk6PyQt4 as v
        else:
            import trussViewerVtk7PyQt5 as v
        print 'Initialize MeshViewer...'
        gui = v.MeshViewer(truss, algo) 
        print 'MeshViewer initialized!'
        gui.start()
    else:
        algo.run() # Runs the algorithm. If you are using the GUI, the GUI will take care of that
        
    # Analytical solution:
    counter = 0
    lower_bound = 0
    upper_bound = 2.2
    step        = 0.01
    u_b = np.zeros((upper_bound-lower_bound)/step)
    for n in np.arange(lower_bound,upper_bound,step):
        u_b[counter] = n
        counter += 1
    P = (E*A*b*b*b)/(2*l0*l0*l0)* np.multiply(u_b,np.multiply(u_b-1,u_b-2))
    u1 = (1+math.sqrt(3)/3)*b;
    u2 = (1-math.sqrt(3)/3)*b;
    Pcrit_analytical = E*A/(2*l0**3)*(u2**3-3*b*u2**2+2*b**2*u2)
    f = open('Node_2_DISPLACEMENTS.ascii', 'r')
    uy2 = []
    for line in f:
        line = line.strip()
        columns = line.split()
        uy2.append(float(columns[1]))
        
    f.close()
    uy2 = np.array(uy2)    
    
    f = open('Lambda.ascii')
    lambda_ = []
    for line in f:
        line = line.strip()
        columns = line.split()
        lambda_.append(float(columns[0]))
        
    f.close()
    lambda_ = np.array(lambda_)
    X = -uy2/b
    Y = 1.5*Pcrit_analytical*lambda_
        
    plot_num = plt.plot(X,Y,'k-.',linewidth=3.0,label='numerical')
    plot_ana = plt.plot(u_b,P,'r-.',linewidth=3.0,label='analytical')
    plt.legend()
    plt.xlabel('displ',fontsize=22)
    plt.ylabel('load P',fontsize=22)
    plt.title('Incremental method',fontsize=26)
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(nogui)