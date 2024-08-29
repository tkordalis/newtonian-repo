#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.5.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/tkordalis/Desktop/SalomeMeshes')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
from math  import *
from numpy import *
import SALOMEDS
from salome.geom import geomtools


# from Geometry_Mesh_Parameters import Rinf,    RSphere1, RSphere2, Sphere_position,\
# 																								 dR_ref1, dR_ref2, \
# R_refinement1_Sphere1, R_refinement1_Sphere2, R_refinement2_Sphere1, R_refinement2_Sphere2, \
# ellipse_position, ellipse_Minor_Radius, ellipse_Major_Radius, outer_ellipse_Major_Radius, outer_ellipse_Minor_Radius, h_s, \
# Main_maxSize_element, Main_minSize_element, Element_size_on_Sphere, Netgen_Params, NumSegmentsOnSphere

from Geometry_Mesh_Parameters import Radius_tank, Height_tank, RSphere1, RSphere2, Sphere_position, dR_ref1, dR_ref2, \
R_refinement1_Sphere1, R_refinement1_Sphere2, R_refinement2_Sphere1, R_refinement2_Sphere2, \
ellipse_position, ellipse_Minor_Radius, ellipse_Major_Radius, outer_ellipse_Major_Radius, outer_ellipse_Minor_Radius, h_s, \
Main_maxSize_element, Main_minSize_element, Element_size_on_Sphere, Netgen_Params, NumSegmentsOnSphere


geompy = geomBuilder.New()

segment_length_from_NumberOfSegments = 3.14159265359/NumSegmentsOnSphere
cb1NumSegments = int(dR_ref1/segment_length_from_NumberOfSegments)


# Base Vectors
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(7, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 7, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 7)




Domain_cut = geompy.MakeCutList( geompy.MakeFaceHW(Height_tank, 2*Radius_tank, 1), \
								[ geompy.MakeTranslation( geompy.MakeDiskR(RSphere1, 1), 0, 0, 0 ), \
								geompy.MakeTranslation( geompy.MakeFaceHW(Height_tank, 2*Radius_tank, 1), 0, -Radius_tank, 0 )], True )


# -------------------------------------------------------
# Domain = geompy.MakeCutList(Cut_1, [Translation_3], True)
# -------------------------------------------------------



Disk_refinement1_Sphere1 = geompy.MakeDiskR(R_refinement1_Sphere1,1)
Disk_refinement2_Sphere1 = geompy.MakeDiskR(R_refinement2_Sphere1,1)

Disk_refinement1_Sphere1_minusSpherePosition = geompy.MakeTranslation(Disk_refinement1_Sphere1, 0, 0, 0)
Disk_refinement2_Sphere1_minusSpherePosition = geompy.MakeTranslation(Disk_refinement2_Sphere1, 0, 0, 0)
[Wire_1] = geompy.ExtractShapes(Disk_refinement1_Sphere1_minusSpherePosition, geompy.ShapeType["WIRE"], True)
[Wire_2] = geompy.ExtractShapes(Disk_refinement2_Sphere1_minusSpherePosition, geompy.ShapeType["WIRE"], True)


vertex_ellipse_position  = geompy.MakeVertex(ellipse_position, 0, 0)

Ellipse_1     = geompy.MakeEllipse(vertex_ellipse_position, None, ellipse_Major_Radius, ellipse_Minor_Radius)
Ellipse_outer = geompy.MakeEllipse(vertex_ellipse_position, None, outer_ellipse_Major_Radius, outer_ellipse_Minor_Radius)



ellipse_position_left  = ellipse_Major_Radius - ellipse_position
ellipse_position_right = ellipse_Major_Radius + ellipse_position

outer_ellipse_position_left  = ellipse_position_left	+ 1.0
outer_ellipse_position_right = ellipse_position_right	+ 1.0



PartitionTool = geompy.MakeFuseList([Wire_1, Wire_2, Ellipse_1, Ellipse_outer], True, True)
Partition_1   = geompy.MakePartition([Domain_cut], [PartitionTool], [], [], geompy.ShapeType["FACE"], 0, [], 0)


def returnPointsFromRotatedAxis( theta_degrees, x_tilt, x_o ):
	x_global 		= []
	theta_radians 	= 3.1415926*(theta_degrees/180.0)

	x_global.append(x_o[0] + cos(theta_radians)*x_tilt[0] + (-sin(theta_radians))*x_tilt[1])
	x_global.append(x_o[1] + sin(theta_radians)*x_tilt[0] +   cos(theta_radians) *x_tilt[1])
	x_global.append(0.0)

	return x_global

def returnIDofShape( theta_degrees, x_tilt, x_o, TypeofShape ):
	idShape     		 =  0.
	Shape_Point_Position = 	returnPointsFromRotatedAxis( theta_degrees, x_tilt, x_o )
	Shape_Point 		 = 	geompy.MakeVertex(*Shape_Point_Position)
	sLineShape  		 = 	geompy.GetShapesNearPoint(Partition_1, Shape_Point, geompy.ShapeType[ TypeofShape ])
	idShape 			 =  geompy.GetSubShapeID(Partition_1,	sLineShape  )
	# a1 = Shape_Point
	# geompy.addToStudy(a1, "a1")

	return idShape

def getgroupSymmetry( mainGroup, refZone, leftOrRight, gid ):
	groupSymmetryDict = {
	"mainGroup"	  : mainGroup,
	"refZone"	  : refZone,
	"leftOrRight" : leftOrRight,
	"id"		  :	gid
	}
	return groupSymmetryDict

def getgroupFace( mainGroup, refZone, gid, obj, union, mesh_obj):
	groupSphereDict = {
	"mainGroup"	  : mainGroup,
	"refZone"	    : refZone,
	"id"		    	: gid,
	"obj"		    	: obj,
	"union"		    : union,
	"mesh_obj" 		: mesh_obj
	}
	return groupSphereDict




# h_s = 0.0001

# # ========================================================================
# # tankWall Group Creation
idtankWall = []
tankWall       = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idtankWall.append(returnIDofShape(0, [ 0, Radius_tank ], [0, 0], "EDGE"))
idtankWall.append( returnIDofShape(0, [ -0.5*Height_tank, +h_s ], [0, 0], "EDGE") ) 
tankWall_union = geompy.UnionIDs(tankWall , idtankWall  )
# geompy.addToStudyInFather(Partition_1, tankWall, "tankWall")



idrightPlane = [] 
rightPlane   = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idrightPlane.append( returnIDofShape(0, [ 0.5*Height_tank, +h_s ], [0, 0], "EDGE") ) 
rightPlane_union = geompy.UnionIDs( rightPlane , idrightPlane )



# idleftPlane = [] 
# leftPlane   = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
# idleftPlane.append( returnIDofShape(0, [ -0.5*Height_tank, +h_s ], [0, 0], "EDGE") ) 
# leftPlane_union = geompy.UnionIDs( leftPlane , idleftPlane )

# # ========================================================================
idSymmetry = []

Symmetry_groups = []


Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb1'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb1'     , 'right', None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb2'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb2'     , 'right', None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'ellipse1', 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'ellipse1', 'right' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'ellipse2', 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'ellipse2', 'right' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'out'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'out'     , 'right' , None) )



origin_of_axes = []
x_tilt 		   = []
theta_degrees  = []

point_list = []

for i in range(len(Symmetry_groups)):

	if Symmetry_groups[i]["mainGroup"] == "Sphere_1_Side":
		origin_of_axes = [0,0]

		if Symmetry_groups[i]["refZone"] == "cb1":
			x_tilt =  [RSphere1+h_s, 0] 

		elif Symmetry_groups[i]["refZone"] == "cb2":
			x_tilt = [R_refinement1_Sphere1+h_s, 0]

		elif Symmetry_groups[i]["refZone"] == "ellipse1":
			x_tilt = [R_refinement2_Sphere1+h_s, 0]

		elif Symmetry_groups[i]["refZone"] == "ellipse2":
			x_tilt = [ellipse_position_left-0+h_s, 0] 

		elif Symmetry_groups[i]["refZone"] == "out":
			x_tilt = [outer_ellipse_position_left-0+h_s, 0] 

	

	if Symmetry_groups[i]["leftOrRight"] == "left":
		theta_degrees = 180
	elif Symmetry_groups[i]["leftOrRight"] == "right":
		theta_degrees = 0

	Symmetry_groups[i]["id"] = returnIDofShape( theta_degrees, x_tilt, origin_of_axes, "EDGE" )
	idSymmetry.append( Symmetry_groups[i]["id"] )

	# print( Symmetry_groups[i] )

	if ( Symmetry_groups[i]["refZone"] == "out" and Symmetry_groups[i]["leftOrRight"] == "left" ):
		SymmetryB1out = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
		SymmetryB1out_union = geompy.UnionIDs( SymmetryB1out, [ Symmetry_groups[i]["id"] ] )
	elif ( Symmetry_groups[i]["refZone"] == "out" and Symmetry_groups[i]["leftOrRight"] == "right" ) :
		SymmetryB2out = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
		SymmetryB2out_union = geompy.UnionIDs( SymmetryB2out, [ Symmetry_groups[i]["id"] ] )





Symmetry   = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Symmetry_union = geompy.UnionIDs(Symmetry , idSymmetry )
# geompy.addToStudyInFather(Partition_1, Symmetry, "Symmetry")
idhorizontalsSphere1 = []
for item in Symmetry_groups:
	if item["refZone"] == "cb1":
			idhorizontalsSphere1.append(item["id"])



horizontalsSphere1       = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
horizontalsSphere1_union = geompy.UnionIDs( horizontalsSphere1, idhorizontalsSphere1 )

geompy.addToStudyInFather(Partition_1, horizontalsSphere1, "horizontalsSphere1")

idSphere1 	  = [returnIDofShape( 1, [RSphere1, 0], [-0,0], "EDGE" )]
Sphere1 	  = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Sphere1_union = geompy.UnionIDs( Sphere1 , idSphere1 )

# print(idSphere1)



idSphere1Ref1     = [returnIDofShape( 1, [RSphere1+dR_ref1, 0], [-0,0], "EDGE" )]
Sphere1Ref1 	    = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Sphere1_unionRef1 = geompy.UnionIDs( Sphere1Ref1 , idSphere1Ref1 )

geompy.addToStudyInFather(Partition_1,Sphere1Ref1,"Sphere1Ref1")
# print(idSphere1Ref1,idSphere2Ref1)

Groups_faces = []

Groups_faces.append( getgroupFace("Sphere1", "Ref_1", None, None, None, None) )
Groups_faces.append( getgroupFace("Sphere1", "Ref_2", None, None, None, None) )
Groups_faces.append( getgroupFace("ellipse", "1"    , None, None, None, None) )
Groups_faces.append( getgroupFace("ellipse", "2"    , None, None, None, None) )

for i in range(len(Groups_faces)):

	if Groups_faces[i]["mainGroup"] == "Sphere1":
		origin_of_axes = [0,0]

		if Groups_faces[i]["refZone"] == "Ref_1":
			x_tilt = [0,RSphere1 + h_s]
		elif Groups_faces[i]["refZone"] == "Ref_2":
			x_tilt = [0,R_refinement1_Sphere1 + h_s]
	elif Groups_faces[i]["mainGroup"] == "ellipse":
		if Groups_faces[i]["refZone"] == "1":
			x_tilt = [0,R_refinement2_Sphere1 + h_s]
		elif Groups_faces[i]["refZone"] == "2":
			x_tilt = [0,ellipse_Minor_Radius + h_s]

	theta_degrees = 0

	Groups_faces[i]["id"] = returnIDofShape( theta_degrees, x_tilt, origin_of_axes, "FACE")

i=0
for item in Groups_faces:
	i = i+1
	item["obj"] = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
	item["union"] = geompy.UnionIDs(item["obj"] , [ item["id" ]] )
	# print(item["mainGroup"], item["refZone"])
	# geompy.addToStudyInFather(Partition_1, item["obj"], str(i))




# # --------------------------- End of Geometry --------------------------- #

# ###
# ### SMESH component
# ###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()


def MeshParameters( Properties, maxsize: float, minsize: float, growthRate: float): 

    params = Properties.Parameters()
    params.SetMaxSize            (  maxsize  )
    params.SetMinSize            (  minsize  )
    params.SetSecondOrder        (   0 )
    params.SetOptimize           (   1 )
    params.SetFineness           (   2 )
    params.SetGrowthRate         (  growthRate )
    params.SetChordalError       (  -1 )
    params.SetChordalErrorEnabled(   0 )
    params.SetUseSurfaceCurvature(   1 )
    params.SetFuseEdges          (   1 )
    params.SetWorstElemMeasure   (   0 )
    params.SetUseDelauney        ( 108 )
    params.SetQuadAllowed        (   0 )
    params.SetCheckChartBoundary (   0 )




Mesh_1 = smesh.Mesh(Partition_1)
# Mesh_1.Segment(geom=Sphere1).StartEndLength (Element_size_on_Sphere,Element_size_on_Sphere,[])
# Mesh_1.Segment(geom=Sphere2).StartEndLength (Element_size_on_Sphere,Element_size_on_Sphere,[])

Mesh_1.Segment(geom=Sphere1).NumberOfSegments(NumSegmentsOnSphere)

Mesh_1.Segment(geom=Sphere1Ref1).NumberOfSegments(NumSegmentsOnSphere)
Mesh_1.Segment(geom=horizontalsSphere1).NumberOfSegments(cb1NumSegments)







NETGEN_1D_2D   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
MeshParameters(NETGEN_1D_2D  ,Main_maxSize_element, Main_minSize_element, 0.1)

Params = []

Groups_faces[0]["mesh_obj"] = Mesh_1.Quadrangle(geom = Groups_faces[0]["obj"]) 

for i in range(len(Groups_faces)-1):
	Params = Netgen_Params[i+1]
	Groups_faces[i+1]["mesh_obj"] = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = Groups_faces[i+1]["obj"])
	MeshParameters(Groups_faces[i+1]["mesh_obj"], Params[0], Params[1], Params[2])

SymmOut_params = Netgen_Params[-1]
Mesh_1.Segment(geom=SymmetryB1out).StartEndLength ( Main_maxSize_element, SymmOut_params[0] )
Mesh_1.Segment(geom=SymmetryB2out).StartEndLength ( SymmOut_params[0], Main_maxSize_element )

tankWall_1	 	=  Mesh_1.GroupOnGeom( tankWall 	,'tankWall'   ,SMESH.EDGE )
Symmetry_1	 	=  Mesh_1.GroupOnGeom( Symmetry 	,'Symmetry'	  ,SMESH.EDGE )
Sphere1_1	 		=  Mesh_1.GroupOnGeom( Sphere1 		,'Bubble1' 	  ,SMESH.EDGE )
rightPlane_1 	=  Mesh_1.GroupOnGeom( rightPlane ,'Ambient' 		,SMESH.EDGE )


Priority_list = []

for item in Groups_faces:
	SubMesh = item["mesh_obj"]
	Priority_list.append(SubMesh.GetSubMesh())

# print(Priority_list)
isDone = Mesh_1.SetMeshOrder( [Priority_list] )

isDone = Mesh_1.Compute()

# isDone = Mesh_1.QuadTo4Tri( )
isDone = Mesh_1.SplitQuadObject( Mesh_1, 1 )



try:
  Mesh_1.ExportUNV( r'./UnBounded.unv' )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')



if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()


import os

try:
    from killSalomeWithPort import killMyPort
    killMyPort(os.getenv('NSPORT'))
except:
    pass

