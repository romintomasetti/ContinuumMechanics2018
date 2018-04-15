# -*- coding: latin-1; -*-
# $Id$

# --- GUI for Trusses to use with VTK 6 and PyQt4 --- #

import sys
import vtk
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

class MeshViewer(QtGui.QWidget):

    app = QtGui.QApplication(sys.argv)
    
    """
    Qt GUI for visu. the output
    """
    def __init__(self, truss, algo): 
        QtGui.QWidget.__init__(self) 
        
        self.truss = truss
        self.algo = algo
        self.algo.dhook = self
        
        self.running = 'init'
        
        print "starting MeshViewer init..."
        
        self.__setupGUI()
        self.__setupVTK()
        
        self.app.lastWindowClosed.connect(self.app.quit)
        self.show() 
        print "MeshViewer ready." 
        
    def closeEvent(self, event):
        self.algo.stopit=True 
        self.running='running' # sort de "while self.running=='pause'"      
        print "GUI killed!"
        QtGui.QWidget.closeEvent(self,event)


    def start(self):
        self.app.exec_()

    def __setupGUI(self):

        self.setWindowTitle("MeshViewer")
        self.resize(800, 600)
        
        # vtk window
        
        self.vtkwidget = QVTKRenderWindowInteractor(self) # "self" sinon, rien ne s'affiche
        self.vtkwidget.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, 
                                                  QtGui.QSizePolicy.Expanding))
        self.vtkwidget.setMinimumSize(QtCore.QSize(300, 300));
        self.vtkwidget.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        
                      
        self.vtkwidget.Initialize() 
                       
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1.0, 1.0, 1.0)     
        self.vtkwidget.GetRenderWindow().AddRenderer(self.renderer) 
        
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.vtkwidget.SetInteractorStyle(style) 


        hbox = QtGui.QHBoxLayout()
        self.setLayout(hbox)
        hbox.addWidget(self.vtkwidget)
        
        pan = QtGui.QFrame()
        pan.setMaximumSize(QtCore.QSize(200,999999))
        hbox.addWidget(pan)
        
        vbox = QtGui.QVBoxLayout()
        pan.setLayout(vbox)
        
        self.startBut = QtGui.QPushButton(self.tr("start!"))
        self.connect(self.startBut, QtCore.SIGNAL("clicked()"), self.startSlot)
        vbox.addWidget(self.startBut)
        
        groupBox = QtGui.QGroupBox("Infos")
        self.steplabel = QtGui.QLabel("step # 0")
        self.loadlabel = QtGui.QLabel("lambda = %2.8f" % 0)
        gbox = QtGui.QVBoxLayout()
        groupBox.setLayout(gbox) 
        gbox.addWidget(self.steplabel)   
        gbox.addWidget(self.loadlabel)   
        vbox.addWidget(groupBox)
        
        vbox.addStretch(1)
        
    def startSlot(self):
        if self.running=='init':  
            self.startBut.setText('Pause') # on demarre et on affiche "pause"
            self.running='running'
            self.algo.run()
            self.startBut.setText("Quit")
            self.running='quit'
        elif self.running=='running':  # on stoppe et on affiche 'continue"
            self.running='pause'
            self.startBut.setText("Continue")
            while self.running=='pause':
                self.app.processEvents(QtCore.QEventLoop.WaitForMoreEvents)
        elif self.running=='pause':
            self.running='running'
            self.startBut.setText("Pause")
        elif self.running=='quit':
            self.app.quit()
        
    def disableStart(self):
        self.startBut.setDisabled(True)
        
    def __setupVTK(self):
        # polydata
        self.__createPolyData()
        self.poly = PolyData(self.polydata)
        self.renderer.AddActor(self.poly.actor)
        self.renderer.AddActor2D(self.poly.pointLabels)
        
        self.resetCamera() 

    def resetCamera(self):
        self.renderer.ResetCamera()
        cam1 = self.renderer.GetActiveCamera()
        # 3D
        if 0:
            cam1.Elevation(35)
            cam1.SetViewUp(0, 1, 0)
            cam1.Azimuth(30)
        #2D
        else:
            cam1.Elevation(0)
            cam1.SetViewUp(0, 1, 0)
            cam1.Azimuth(0)
            self.renderer.ResetCameraClippingRange()        
        
    def display(self, step, lamda):
        
        self.steplabel.setText("step # %d" % step)
        self.loadlabel.setText("lambda = %2.8f" % lamda)
        
        self.points.Reset()
        self.vertices.Reset()
        self.lines.Reset()
        self.scalars.Reset()
        
        # points
        nmap={}
        i=0
        for nod in self.truss.nodes:
            nmap[nod.nb]=i
            self.points.InsertPoint(i, nod.x, nod.y, 0.0)
            self.scalars.InsertNextValue(nod.nb)
            # node cells
            vertex = vtk.vtkVertex()
            vertex.GetPointIds().SetId(0, i)
            self.vertices.InsertNextCell(vertex)
            i+=1
        self.points.Modified()
        self.vertices.Modified()
        
        for bar in self.truss.bars:
            line = vtk.vtkLine()
            ids = line.GetPointIds()
            ids.SetNumberOfIds(2)
            for j in range(len(bar.nodes)):
                ids.SetId(j, nmap[bar.nodes[j].nb])
            self.lines.InsertNextCell(line) 
        self.lines.Modified()
        
        self.polydata.Modified()
        
        self.render()
    
    def ragequit(self):
        print "rage quit!"
        self.algo.stopit=True
        self.app.quit()
        
    def render(self): 
        # draw the scene
        self.vtkwidget.Render()
        self.app.processEvents()
    
    def refresh(self):
        self.app.processEvents()
        
    def __createPolyData(self):
        print 'creating vtkPolyData...'
        self.points = vtk.vtkPoints()
        self.polydata  = vtk.vtkPolyData()
        self.vertices = vtk.vtkCellArray()
        self.lines = vtk.vtkCellArray()
        
        self.polydata.SetPoints(self.points)
        self.polydata.SetVerts(self.vertices)
        self.polydata.SetLines(self.lines)
        
        # points
        self.scalars = vtk.vtkFloatArray()
        self.scalars.SetNumberOfComponents(1)
        self.polydata.GetPointData().SetScalars(self.scalars)
     
        nmap={}
        i=0
        for nod in self.truss.nodes:
            nmap[nod.nb]=i
            self.points.InsertPoint(i, nod.x, nod.y, 0.0)
            self.scalars.InsertNextValue(nod.nb)
            # node cells
            vertex = vtk.vtkVertex()
            vertex.GetPointIds().SetId(0, i)
            self.vertices.InsertNextCell(vertex)
            i+=1
        self.points.Modified()
        self.vertices.Modified()
        
        for bar in self.truss.bars:
            line = vtk.vtkLine()
            ids = line.GetPointIds()
            ids.SetNumberOfIds(2)
            for j in range(len(bar.nodes)):
                ids.SetId(j, nmap[bar.nodes[j].nb])
            self.lines.InsertNextCell(line) 
        self.lines.Modified()

class PolyData:
    def __init__(self, polydata):
        
        self.polydata = polydata
        
        self.mapper = vtk.vtkPolyDataMapper() 
        self.mapper.ImmediateModeRenderingOff()
        self.mapper.ScalarVisibilityOff()
        self.mapper.SetInputDataObject(polydata)
        self.actor = vtk.vtkActor()
            
        self.actor.GetProperty().SetPointSize(7.0)
        self.actor.GetProperty().SetColor(0.,0.,0.)
        self.actor.SetMapper(self.mapper)
        
        self.labelmapper = vtk.vtkLabeledDataMapper()
        self.labelmapper.SetLabelFormat(" %g")
        self.labelmapper.SetLabelModeToLabelScalars()
        self.labelmapper.SetInputDataObject(polydata)
        prp = self.labelmapper.GetLabelTextProperty()
        prp.SetFontSize(16)
        prp.SetItalic(False)
        prp.SetBold(True)
        prp.SetColor(1.0,0.0,0.0)
        prp.SetLineOffset(-2)
        
        self.pointLabels = vtk.vtkActor2D()
        self.pointLabels.SetMapper(self.labelmapper)