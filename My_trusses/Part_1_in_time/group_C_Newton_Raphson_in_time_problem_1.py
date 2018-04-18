import sys
sys.path.insert(0, '../')
from node import *
from bar import *
from truss import *
from nonLinAlgo import *

from optparse import OptionParser
import vtk
import math

def main(nogui):
    
    #Geometry
    a = 2.0
    b = 1.0
    #sys.exit()
    l0 = math.sqrt(a**2 + b**2)
    tol = 1e-9
    nItMax = 1000


    
    #Nodes
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, a, b)
    
    #Material properties
    E = 70e9
    A = 0.01
    
    #Critical load
    qcr = (math.sqrt(3)*E*A*(b**3))/(9.0*(l0**3)) #It's the critical load divided by 2 since we consider just one bar!
    #We go a little bit higher than the critical load:
    F_max = -qcr*1.2
    print F_max
    #Max time before stopping:
    t_max = 2.0
    #Time step:
    dt = 0.005
    #Period:
    T = 1
    
    if(dt > T or dt > t_max or T > t_max):
        print "Please check your parameters"
        sys.exit()
    
    #Bars
    bar1 = Bar(1, [node1, node2], E, A)
    
    #Truss
    truss = Truss()
    truss.addNode(node1)
    truss.addNode(node2)
    truss.addBar(bar1)
    truss.setNodeRows() # ATTENTION: this line is mandatory at the end of the definition of the truss!
    
    #BCs
    truss.fix(node1,'x')
    truss.fix(node1,'y')
    truss.fix(node2,'x')
    
    #Non-linear algorithm
    algo = NewtonRaphsonAlgorithm_inTime(truss, tol, nItMax,dt, t_max, F_max, T)
    algo.run()

if __name__ == '__main__':
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(nogui)