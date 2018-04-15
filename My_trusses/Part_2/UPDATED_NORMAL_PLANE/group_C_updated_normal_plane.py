import sys
sys.path.insert(0, '../')
from node import *
from bar import *
from truss import *
from nonLinAlgo import *

import matplotlib.pyplot as plt

from optparse import OptionParser
import vtk
import math

def main(nogui):
    
    # Geometry:
    a     = 0.75
    b     = 0.25
    alpha = 30*math.pi/180
    
    # Material:
    E_1 = 70e9
    E_2 = 70e9
    A_1 = 25e-4
    A_2 = 25e-4
    
    # Misc.
    l0 = math.sqrt(a**2 + b**2)
    toll = 1e-6
    nItMax = 60
    
    #Nodes
    node1 = Node(1 ,0.0, 0.0)
    node2 = Node(2 ,a*math.cos(alpha)                     , a*math.sin(alpha))
    node3 = Node(3, a*math.cos(alpha)-b                   , a*math.sin(alpha))
    node4 = Node(4, 2*a*math.cos(alpha)-b                   , 0.0         )
    
    #Bars
    bar1 = Bar(1, [node1, node2], E_1, A_1)
    bar2 = Bar(2, [node2, node3], E_2, A_2)
    bar3 = Bar(3, [node3, node4], E_1, A_1)
    
    #Truss
    truss = Truss()
    truss.addNode(node1)
    truss.addNode(node2)
    truss.addNode(node3)
    truss.addNode(node4)
    truss.addBar(bar1)
    truss.addBar(bar2)
    truss.addBar(bar3)
    truss.setNodeRows() # ATTENTION: this line is mandatory at the end of the definition of the truss!
    
    #BCs
    truss.fix(node1,'x')
    truss.fix(node1,'y')
    truss.fix(node4,'x')
    truss.fix(node4,'y')
    
    #Critical load
    qcr = 25*(math.sqrt(3)*E_2*A_2*(b**3))/(9.0*(l0**3)) #It's the critical load divided by 2 since we consider just one bar!

    #Loads
    node2.applyLoad('y', -1.5*qcr)
    node3.applyLoad('y', -1.5*qcr)
    
    #Use a copy of the truss.
    truss_test = truss;
    
    #Non-linear algorithm
    dlambda_0 = 1e-3
    Id_0 = 10
    psi_0 = 1e-10
    while dlambda_0 <= 1e-2:
        Id = Id_0
        while Id <= 1000:
            psi = psi_0
            while psi <= 1:
                dlamda0 = dlambda_0
                algo = UpdatedNormalPlaneArcLengthAlgorithm(truss_test, toll, nItMax,dlamda0, psi, Id)
                #algo = NewtonRaphsonAlgorithm(truss,toll,nItMax,dlamda0)
                #algo = ArcLengthAlgorithm(truss,toll,nItMax,dlamda0,psi,Id)
                #algo = IncrementalAlgorithm(truss,dlamda0)
                print "Starting with (",dlamda0,",",Id,",",psi,")."
                returned_value = algo.run()
                if not truss.nodes[0].x == 0:
                    sys.exit("The algo modified the initial truss!")
                if returned_value == -7:
                    print "It failed with (",dlamda0,",",Id,",",psi,")."
                else:
                    print "Success with (",dlamda0,",",Id,",",psi,")."
                    sys.exit("SUCCESS")
                psi *= 10
                sys.exit()
                #End of loop on psi
            Id += 10
            #End of loop on Id
        dlambda_0 *= 10
        #End of loop on dlambda
            

if __name__ == '__main__':
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(nogui)