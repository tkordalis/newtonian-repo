#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.9.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/tkordalis/Documents/Salome')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

Rbubble1 = 1.0

dR_ref1 = 0.050
dR_ref2 = 0.150
R_refinement1_Bubble1 = Rbubble1 + dR_ref1
R_refinement2_Bubble1 = Rbubble1 + dR_ref2


O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Disk_1 = geompy.MakeDiskR(500, 1)
Disk_2 = geompy.MakeDiskR(Rbubble1, 1)
Cut_1 = geompy.MakeCutList(Disk_1, [Disk_2], True)
Disk_3 = geompy.MakeDiskR(R_refinement1_Bubble1, 1)
[Wire_1] = geompy.ExtractShapes(Disk_3, geompy.ShapeType["WIRE"], True)
Disk_4 = geompy.MakeDiskR(R_refinement2_Bubble1, 1)
[Wire_2] = geompy.ExtractShapes(Disk_4, geompy.ShapeType["WIRE"], True)
Disk_5 = geompy.MakeDiskR(3.1, 1)
[Wire_3] = geompy.ExtractShapes(Disk_5, geompy.ShapeType["WIRE"], True)
Disk_6 = geompy.MakeDiskR(4.1, 1)
[Wire_4] = geompy.ExtractShapes(Disk_6, geompy.ShapeType["WIRE"], True)
Face_1 = geompy.MakeFaceHW(1100, 1100, 1)
Translation_1 = geompy.MakeTranslation(Face_1, 0, -550, 0)
Partition_2 = geompy.MakePartition([Cut_1], [Wire_2, Wire_3, Wire_4, Wire_1], [], [], geompy.ShapeType["FACE"], 0, [], 0)
Cut_2 = geompy.MakeCutList(Partition_2, [Translation_1], True)
Outflow = geompy.CreateGroup(Cut_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Outflow, [5])
Symmetry = geompy.CreateGroup(Cut_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Symmetry, [44, 30, 35, 23, 19, 40, 36, 31, 8, 12])
Bubble1 = geompy.CreateGroup(Cut_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bubble1, [42])
Ref1 = geompy.CreateGroup(Cut_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Ref1, [38])
Ref2 = geompy.CreateGroup(Cut_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Ref2, [14])
Ref3 = geompy.CreateGroup(Cut_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Ref3, [25])
Ref4 = geompy.CreateGroup(Cut_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Ref4, [33])
[Outflow, Symmetry, Bubble1, Ref1, Ref2, Ref3, Ref4] = geompy.GetExistingSubObjects(Cut_2, False)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Disk_1, 'Disk_1' )
geompy.addToStudy( Disk_2, 'Disk_2' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudyInFather( Disk_4, Wire_2, 'Wire_2' )
geompy.addToStudyInFather( Disk_5, Wire_3, 'Wire_3' )
geompy.addToStudyInFather( Disk_6, Wire_4, 'Wire_4' )
geompy.addToStudyInFather( Disk_3, Wire_1, 'Wire_1' )
geompy.addToStudy( Partition_2, 'Partition_2' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_2, 'Cut_2' )
geompy.addToStudyInFather( Cut_2, Outflow, 'Outflow' )
geompy.addToStudyInFather( Cut_2, Symmetry, 'Symmetry' )
geompy.addToStudy( Disk_4, 'Disk_4' )
geompy.addToStudy( Disk_5, 'Disk_5' )
geompy.addToStudy( Disk_6, 'Disk_6' )
geompy.addToStudyInFather( Cut_2, Bubble1, 'Bubble1' )
geompy.addToStudyInFather( Cut_2, Ref1, 'Ref1' )
geompy.addToStudyInFather( Cut_2, Ref2, 'Ref2' )
geompy.addToStudy( Disk_3, 'Disk_3' )
geompy.addToStudyInFather( Cut_2, Ref3, 'Ref3' )
geompy.addToStudyInFather( Cut_2, Ref4, 'Ref4' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Cut_2,'Mesh_1')
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 30 )
NETGEN_2D_Parameters_1.SetMinSize( 0.1 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 21852 )
NETGEN_2D_Parameters_1.SetUseDelauney( 224 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 3 )

Element_size_on_Bubble = 0.002

Netgen_Params1 = []
Netgen_Params2 = []
Netgen_Params3 = []
Netgen_Params4 = []

Netgen_Params1 = [ 4*Element_size_on_Bubble,  4*Element_size_on_Bubble, 0.5]
Netgen_Params2 = [ 8*Element_size_on_Bubble,  8*Element_size_on_Bubble, 0.6]
Netgen_Params3 = [      0.05       ,        0.05      , 0.5]
Netgen_Params4 = [     0.12       ,        0.1       , 0.1]

Regular_1D = Mesh_1.Segment(geom=Bubble1)
# Start_and_End_Length_1 = Regular_1D.StartEndLength(0.006,0.006,[])
Start_and_End_Length_1 = Regular_1D.StartEndLength(Element_size_on_Bubble,Element_size_on_Bubble,[])

Start_and_End_Length_1.SetObjectEntry( "0:1:1:18" )
NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Ref1)
Sub_mesh_2 = NETGEN_1D_2D_1.GetSubMesh()
NETGEN_2D_Parameters_2 = NETGEN_1D_2D_1.Parameters()

# NETGEN_2D_Parameters_2.SetMaxSize( 0.012 )
# NETGEN_2D_Parameters_2.SetMinSize( 0.012 )

NETGEN_2D_Parameters_2.SetMaxSize( Netgen_Params1[0] )
NETGEN_2D_Parameters_2.SetMinSize( Netgen_Params1[1] )

NETGEN_2D_Parameters_2.SetSecondOrder( 0 )
NETGEN_2D_Parameters_2.SetOptimize( 1 )
NETGEN_2D_Parameters_2.SetFineness( 5 )
NETGEN_2D_Parameters_2.SetGrowthRate( Netgen_Params1[2])
NETGEN_2D_Parameters_2.SetNbSegPerEdge( 15 )
NETGEN_2D_Parameters_2.SetNbSegPerRadius( 2 )
NETGEN_2D_Parameters_2.SetChordalError( -1 )
NETGEN_2D_Parameters_2.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_2.SetFuseEdges( 1 )
NETGEN_2D_Parameters_2.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_2.SetWorstElemMeasure( 21852 )
NETGEN_2D_Parameters_2.SetUseDelauney( 224 )
NETGEN_2D_Parameters_2.SetCheckChartBoundary( 3 )
NETGEN_1D_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Ref2)
Sub_mesh_3 = NETGEN_1D_2D_2.GetSubMesh()
NETGEN_2D_Parameters_3 = NETGEN_1D_2D_2.Parameters()

# NETGEN_2D_Parameters_3.SetMaxSize( 0.024 )
# NETGEN_2D_Parameters_3.SetMinSize( 0.024 )

NETGEN_2D_Parameters_3.SetMaxSize( Netgen_Params2[0] )
NETGEN_2D_Parameters_3.SetMinSize( Netgen_Params2[1] )

NETGEN_2D_Parameters_3.SetSecondOrder( 0 )
NETGEN_2D_Parameters_3.SetOptimize( 1 )
NETGEN_2D_Parameters_3.SetFineness( 5 )
NETGEN_2D_Parameters_3.SetGrowthRate( Netgen_Params2[2] )
NETGEN_2D_Parameters_3.SetNbSegPerEdge( 15 )
NETGEN_2D_Parameters_3.SetNbSegPerRadius( 2 )
NETGEN_2D_Parameters_3.SetChordalError( -1 )
NETGEN_2D_Parameters_3.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_3.SetFuseEdges( 1 )
NETGEN_2D_Parameters_3.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_3.SetWorstElemMeasure( 21852 )
NETGEN_2D_Parameters_3.SetUseDelauney( 224 )
NETGEN_2D_Parameters_3.SetCheckChartBoundary( 3 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_2, Sub_mesh_3 ] ])
NETGEN_1D_2D_3 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Ref3)
Sub_mesh_4 = NETGEN_1D_2D_3.GetSubMesh()
NETGEN_2D_Parameters_4 = NETGEN_1D_2D_3.Parameters()

# NETGEN_2D_Parameters_4.SetMaxSize( 0.05 )
# NETGEN_2D_Parameters_4.SetMinSize( 0.05 )

NETGEN_2D_Parameters_4.SetMaxSize( Netgen_Params3[0] )
NETGEN_2D_Parameters_4.SetMinSize( Netgen_Params3[1] )

NETGEN_2D_Parameters_4.SetSecondOrder( 0 )
NETGEN_2D_Parameters_4.SetOptimize( 1 )
NETGEN_2D_Parameters_4.SetFineness( 5 )
NETGEN_2D_Parameters_4.SetGrowthRate( Netgen_Params3[2] )
NETGEN_2D_Parameters_4.SetNbSegPerEdge( 15 )
NETGEN_2D_Parameters_4.SetNbSegPerRadius( 2 )
NETGEN_2D_Parameters_4.SetChordalError( -1 )
NETGEN_2D_Parameters_4.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_4.SetFuseEdges( 1 )
NETGEN_2D_Parameters_4.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_4.SetWorstElemMeasure( 21852 )
NETGEN_2D_Parameters_4.SetUseDelauney( 224 )
NETGEN_2D_Parameters_4.SetCheckChartBoundary( 3 )
isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_2, Sub_mesh_3, Sub_mesh_4 ] ])
NETGEN_1D_2D_4 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=Ref4)
Sub_mesh_5 = NETGEN_1D_2D_4.GetSubMesh()
NETGEN_2D_Parameters_5 = NETGEN_1D_2D_4.Parameters()

# NETGEN_2D_Parameters_5.SetMaxSize( 0.1 )
# NETGEN_2D_Parameters_5.SetMinSize( 0.1 )

NETGEN_2D_Parameters_5.SetMaxSize( Netgen_Params4[0] )
NETGEN_2D_Parameters_5.SetMinSize( Netgen_Params4[1] )

NETGEN_2D_Parameters_5.SetSecondOrder( 0 )
NETGEN_2D_Parameters_5.SetOptimize( 1 )
NETGEN_2D_Parameters_5.SetFineness( 5 )
NETGEN_2D_Parameters_5.SetGrowthRate( Netgen_Params4[2])
NETGEN_2D_Parameters_5.SetNbSegPerEdge( 15 )
NETGEN_2D_Parameters_5.SetNbSegPerRadius( 2 )
NETGEN_2D_Parameters_5.SetChordalError( -1 )
NETGEN_2D_Parameters_5.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_5.SetFuseEdges( 1 )
NETGEN_2D_Parameters_5.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_5.SetWorstElemMeasure( 21852 )
NETGEN_2D_Parameters_5.SetUseDelauney( 224 )
NETGEN_2D_Parameters_5.SetCheckChartBoundary( 3 )
Outflow_1 = Mesh_1.GroupOnGeom(Outflow,'Outflow',SMESH.EDGE)
Symmetry_1 = Mesh_1.GroupOnGeom(Symmetry,'Symmetry',SMESH.EDGE)
Bubble1_1 = Mesh_1.GroupOnGeom(Bubble1,'Bubble1',SMESH.EDGE)

isDone = Mesh_1.SetMeshOrder( [ [ Sub_mesh_2, Sub_mesh_3, Sub_mesh_4, Sub_mesh_5 ] ])
isDone = Mesh_1.Compute()
# [ Outflow_1, Symmetry_1, Bubble1_1 ] = Mesh_1.GetGroups()
try:
  Mesh_1.ExportUNV( r'./UnBounded.unv' )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')
Sub_mesh_1 = Regular_1D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Start_and_End_Length_1, 'Start and End Length_1')
smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 2D Parameters_2')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(NETGEN_2D_Parameters_5, 'NETGEN 2D Parameters_5')
smesh.SetName(NETGEN_2D_Parameters_3, 'NETGEN 2D Parameters_3')
smesh.SetName(Sub_mesh_5, 'Sub-mesh_5')
smesh.SetName(NETGEN_2D_Parameters_4, 'NETGEN 2D Parameters_4')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Outflow_1, 'Outflow')
smesh.SetName(Bubble1_1, 'Bubble1')
smesh.SetName(Symmetry_1, 'Symmetry')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
