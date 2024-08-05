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
Domain_cut = geompy.MakeCutList( geompy.MakeFaceHW(Height_tank, 2*Radius_tank, 1), \
								[ geompy.MakeTranslation( geompy.MakeDiskR(RSphere1, 1), 0, 0, 0 ), \
								geompy.MakeTranslation( geompy.MakeFaceHW(Height_tank, 2*Radius_tank, 1), 0, -Radius_tank, 0 )], True )

geompy.addToStudy(Domain_cut,"domain")

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
	sLineShape  		 = 	geompy.GetShapesNearPoint(Domain_cut, Shape_Point, geompy.ShapeType[ TypeofShape ])
	idShape 			 =  geompy.GetSubShapeID(Domain_cut,	sLineShape  )
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

idtankWall = []
tankWall       = geompy.CreateGroup(Domain_cut, geompy.ShapeType["EDGE"])
idtankWall.append(returnIDofShape(0, [ 0, Radius_tank ], [0, 0], "EDGE"))
idtankWall.append( returnIDofShape(0, [ -0.5*Height_tank, +h_s ], [0, 0], "EDGE") ) 
tankWall_union = geompy.UnionIDs(tankWall , idtankWall  )

idSymmetry = []
Symmetry   = geompy.CreateGroup(Domain_cut, geompy.ShapeType["EDGE"])
idSymmetry.append(returnIDofShape(0, [ -RSphere1-h_s, 0 ], [0, 0], "EDGE"))
idSymmetry.append(returnIDofShape(0, [  RSphere1+h_s, 0 ], [0, 0], "EDGE"))
Symmetry_union = geompy.UnionIDs(Symmetry , idSymmetry )

idSphere1 = []
Sphere1 	  = geompy.CreateGroup(Domain_cut, geompy.ShapeType["EDGE"])
idSphere1 	  = [returnIDofShape( 1, [RSphere1, 0], [-0,0], "EDGE" )]
Sphere1_union = geompy.UnionIDs( Sphere1 , idSphere1 )

idAmbient = []
Ambient 	  = geompy.CreateGroup(Domain_cut, geompy.ShapeType["EDGE"])
idAmbient 	  = [returnIDofShape( 0, [0.5*Height_tank, h_s], [-0,0], "EDGE" )]
Ambient_union = geompy.UnionIDs( Ambient , idAmbient )

geompy.addToStudyInFather(Domain_cut,tankWall,"tankWall")
geompy.addToStudyInFather(Domain_cut,Symmetry,"Symmetry")
geompy.addToStudyInFather(Domain_cut,Sphere1,"Bubble1")
geompy.addToStudyInFather(Domain_cut,Ambient,"Ambient")


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


Mesh_1 = smesh.Mesh(Domain_cut)

Mesh_1.Segment(geom=Sphere1).NumberOfSegments(NumSegmentsOnSphere)
# Mesh_1.Segment(geom=Ambient).NumberOfSegments(5)

NETGEN_1D_2D   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
MeshParameters(NETGEN_1D_2D  ,Main_maxSize_element, Main_minSize_element, 0.2)

tankWall_1	 	=  Mesh_1.GroupOnGeom( tankWall 	,'tankWall'   ,SMESH.EDGE )
Symmetry_1	 	=  Mesh_1.GroupOnGeom( Symmetry 	,'Symmetry'	  ,SMESH.EDGE )
Sphere1_1	 	=  Mesh_1.GroupOnGeom( Sphere1 		,'Bubble1' 	  ,SMESH.EDGE )
rightPlane_1 	=  Mesh_1.GroupOnGeom( Ambient ,'Ambient' 		,SMESH.EDGE )

isDone = Mesh_1.Compute()

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
