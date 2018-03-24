# -*- coding: utf-8 -*-

import vtk
from vtk.util.colors import banana
 
filename = "ATR72_SI_MTOW_Control_1.stl"

# get data

reader = vtk.vtkSTLReader()
reader.SetFileName(filename)

arrow = vtk.vtkArrowSource()
#arrow.SetTipLength(10)
#arrow.SetShaftRadius(2)

# buile mappers
 
mapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(reader.GetOutput())
else:
    mapper.SetInputConnection(reader.GetOutputPort())

arrowMapper = vtk.vtkPolyDataMapper()
arrowMapper.SetInputConnection(arrow.GetOutputPort())

# create the actors

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(banana)
#actor.RotateZ(45.0)

arrowActor=vtk.vtkActor()
arrowActor.SetMapper(arrowMapper)

# add a transform
transform=vtk.vtkTransform()
transform.Translate(0.0,2.0,0.0)
arrowActor.SetUserTransform(transform)
 
# Create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

#arrowRen= vtk.vtkRenderer()
#renWin.AddRenderer(arrowRen)
 
# Create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
 
# Assign actor to the renderer
ren.AddActor(actor)
#arrowRen.AddActor(arrowActor)
ren.AddActor(arrowActor)
renWin.SetSize(500,500)
ren.GetActiveCamera().Zoom(1)
ren.GetActiveCamera().SetFocalDisk(0.5)
renWin.Render()

# try saving the current scene (must launch renWin.Render() first)
writer=vtk.vtkPNGWriter() #define the extension writer

#define a filter
ImFilter=vtk.vtkWindowToImageFilter()
ImFilter.SetInput(renWin)
ImFilter.SetScale(1)
ImFilter.SetInputBufferTypeToRGBA()
#ImFilter.SetInputBufferTypeToRGB() # alternative

#Read the buffer
ImFilter.ReadFrontBufferOff()
ImFilter.Update()

#File name, link and save
writer.SetFileName("essai.png")
writer.SetInputConnection(ImFilter.GetOutputPort())
writer.Write()
 
# Enable user interface interactor
iren.Initialize()

iren.Start()

