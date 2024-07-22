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
from Geometry_Mesh_Parameters import *

# first attempt will be without the wall thickness
# up to now i have seen that the small thickness of the wall creates problems to the elements
# that exist on that small wall part, since the are trapped and become squeezed
# in the second attemp i will ad also this part and examine it more carefully

geompy = geomBuilder.New()



Domain_translated = geompy.MakeTranslation( geompy.MakeFaceHW(cavity_edge, cavity_edge, 1), 0.5*cavity_edge, 0.5*cavity_edge, 0 )
								


# x_tilt of type list [a,b,c]
def returnIDofShape( x_tilt, TypeofShape ):
	idShape     		 =  0.
	Shape_Point 		 = 	geompy.MakeVertex(*x_tilt)
	sLineShape  		 = 	geompy.GetShapesNearPoint(Domain_translated, Shape_Point, geompy.ShapeType[ TypeofShape ])
	idShape 			 =  geompy.GetSubShapeID(Domain_translated,	sLineShape  )
	
	return idShape


idTopWall = []
TopWall = geompy.CreateGroup(Domain_translated, geompy.ShapeType["EDGE"])
idTopWall.append( returnIDofShape( [          h_s    ,     cavity_edge, 0] , "EDGE" ) )
Wall_union = geompy.UnionIDs( TopWall, idTopWall )


idBottomWall = []
BottomWall = geompy.CreateGroup(Domain_translated, geompy.ShapeType["EDGE"])
idBottomWall.append( returnIDofShape( [             cavity_edge,     h_s, 0] , "EDGE" ) )
BottomWall_union = geompy.UnionIDs( BottomWall, idBottomWall )

idRightWall = []
RightWall = geompy.CreateGroup(Domain_translated, geompy.ShapeType["EDGE"])
idRightWall.append( returnIDofShape( [             h_s,     0, 0] , "EDGE" ) )
RightWall_union = geompy.UnionIDs( RightWall, idRightWall )

idLeftWall = []
LeftWall = geompy.CreateGroup(Domain_translated, geompy.ShapeType["EDGE"])
idLeftWall.append( returnIDofShape( [             0,     h_s, 0] , "EDGE" ) )
LeftWall_union = geompy.UnionIDs( LeftWall, idLeftWall )




# # --------------------------- End of Geometry --------------------------- #

# ###
# ### SMESH component
# ###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Domain_translated,'Mesh_1') 
Regular_1D = Mesh_1.Segment()
Number_of_Segments_1 = Regular_1D.NumberOfSegments(NumSegmentsOnCavityEdge)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
isDone = Mesh_1.Compute()

# isDone = Mesh_1.SplitQuadObject( Mesh_1, 2 )
# isDone = Mesh_1.QuadTo4Tri()


BottomWall_1    				=  Mesh_1.GroupOnGeom( BottomWall 	,'BottomWall' 	 ,SMESH.EDGE )
TopWall_1 		                =  Mesh_1.GroupOnGeom( TopWall     ,'TopWall'      ,SMESH.EDGE )
LeftWall_1                      =  Mesh_1.GroupOnGeom( LeftWall     ,'LeftWall'      ,SMESH.EDGE )
RightWall_1                      =  Mesh_1.GroupOnGeom( RightWall     ,'RightWall'      ,SMESH.EDGE )



isDone = Mesh_1.Compute()


try:
  Mesh_1.ExportUNV( r'./UnBounded.unv' )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')



if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()


import os

# try:
#     from killSalomeWithPort import killMyPort
#     killMyPort(os.getenv('NSPORT'))
# except:
#     pass
