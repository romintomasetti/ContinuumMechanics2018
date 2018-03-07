from node import *
from bar import *
from truss import *
from nonLinAlgo import *

from optparse import OptionParser
import vtk

def main(nogui):
    
    #Geometry
    a = 2.0
    b = 1.0
    l0 = math.sqrt(a**2 + b**2)
    
    #Nodes
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, a, b)
    
    #Material properties
    E = 70e9
    A = 0.01
    
    #Critical load
    qcr = (math.sqrt(3)*E*A*(b**3))/(9.0*(l0**3)) #It's the critical load divided by 2 since we consider just one bar!
    
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
    
    #Loads
    node2.applyLoad('y', -1.5*qcr)
    
    #Non-linear algorithm
    dlamda0 = 0.01
    algo = IncrementalAlgorithm(truss, dlamda0)
    
    #GUI --> if you want to run a test without GUI: in a command window type (without the quotation marks) 'python test.py --nogui'
    if not nogui:
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

if __name__ == '__main__':
    
    parser=OptionParser()
    parser.add_option("--nogui", action="store_true",
                        help="Specify if we need to use the GUI", dest="nogui", default=False)

    (options, args)=parser.parse_args()
    
    nogui = options.nogui
    
    main(nogui)