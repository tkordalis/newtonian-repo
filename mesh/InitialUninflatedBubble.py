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



Domain_cut = geompy.MakeCutList( geompy.MakeTranslation( geompy.MakeFaceHW(Height_tank, Radius_tank, 1), 0.5*Height_tank-Height_needle, 0.5*Radius_tank, 0 ), \
								[ geompy.MakeTranslation( geompy.MakeFaceHW(Height_needle, Radius_needle, 1), -0.5*Height_needle, 0.5*Radius_needle, 0 ) ], True )


# geompy.addToStudy(geompy.MakeLineTwoPnt( geompy.MakeVertex(        0,        Radius_needle-Needle_thickness, 0 ) , geompy.MakeVertex( dR_ref, Radius_needle-Needle_thickness, 0) ), 'a')

partition_line_list = []
partition_line_list = [ geompy.MakeLineTwoPnt( geompy.MakeVertex(             0,        Radius_needle, 0 ) , geompy.MakeVertex(      0, Radius_needle + dR_ref, 0) ) ,
						geompy.MakeLineTwoPnt( geompy.MakeVertex(        dR_ref,        Radius_needle, 0 ) , geompy.MakeVertex( dR_ref, Radius_needle + dR_ref, 0) ) ,
						geompy.MakeLineTwoPnt( geompy.MakeVertex(        dR_ref,        Radius_needle, 0 ) , geompy.MakeVertex( dR_ref,                     0 , 0) ) ,
            geompy.MakeLineTwoPnt( geompy.MakeVertex(        0,        Radius_needle-Needle_thickness, 0 ) , geompy.MakeVertex( dR_ref, Radius_needle-Needle_thickness, 0) ),
						geompy.MakeLineTwoPnt( geompy.MakeVertex(-Height_needle, Radius_needle + dR_ref,0) , geompy.MakeVertex(      0, Radius_needle + dR_ref, 0) ) ,
						geompy.MakeLineTwoPnt( geompy.MakeVertex(             0, Radius_needle + dR_ref,0) , geompy.MakeVertex( dR_ref, Radius_needle + dR_ref, 0) ) ,
            geompy.MakeLineTwoPnt( geompy.MakeVertex(             0,          Radius_needle,0) , geompy.MakeVertex( dR_ref,          Radius_needle, 0) ) ,
						geompy.MakeLineTwoPnt( geompy.MakeVertex( -Height_needle, Radius_needle + dR_ref + dH_ref, 0), geompy.MakeVertex(  dR_ref + dH_ref, Radius_needle + dR_ref + dH_ref, 0) ) ,
            geompy.MakeLineTwoPnt( geompy.MakeVertex( dR_ref + dH_ref, Radius_needle + dR_ref + dH_ref, 0) , geompy.MakeVertex( dR_ref + dH_ref, 0, 0) )]

partition_tool = geompy.MakeFuseList( partition_line_list, True, True)


Partition_1 = geompy.MakePartition([Domain_cut], [partition_tool], [], [], geompy.ShapeType["FACE"], 0, [], 0)

geompy.addToStudy(Partition_1, "Partition_1")

Height_tank = Height_tank - Height_needle

# x_tilt of type list [a,b,c]
def returnIDofShape( x_tilt, TypeofShape ):
	idShape     		 =  0.
	Shape_Point 		 = 	geompy.MakeVertex(*x_tilt)
	sLineShape  		 = 	geompy.GetShapesNearPoint(Partition_1, Shape_Point, geompy.ShapeType[ TypeofShape ])
	idShape 			 =  geompy.GetSubShapeID(Partition_1,	sLineShape  )
	
	return idShape


idWall = []
Wall = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idWall.append( returnIDofShape( [              -h_s,     Radius_needle, 0] , "EDGE" ) )
idWall.append( returnIDofShape( [ Height_tank,     h_s, 0] , "EDGE" ) )
idWall.append( returnIDofShape( [ 0,     Radius_needle-Needle_thickness+h_s, 0] , "EDGE" ) )
# idWall.append( returnIDofShape( [    -Height_needle, Radius_needle+h_s, 0] , "EDGE" ) )
# idWall.append( returnIDofShape( [    -Height_needle, Radius_needle+dR_ref+h_s, 0] , "EDGE" ) )
# idWall.append( returnIDofShape( [    -Height_needle,   Radius_tank-h_s, 0] , "EDGE" ) )
idWall.append( returnIDofShape( [-Height_needle+h_s,       Radius_tank, 0] , "EDGE" ) )
Wall_union = geompy.UnionIDs( Wall, idWall )


idSymmetry = []
Symmetry = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])

idSymmetry.append( returnIDofShape( [             h_s,     0, 0] , "EDGE" ) )
idSymmetry.append( returnIDofShape( [ dR_ref + dH_ref - h_s,     0, 0] , "EDGE" ) )
idSymmetry.append( returnIDofShape( [ Height_tank-h_s,     0, 0] , "EDGE" ) )
Symmetry_union = geompy.UnionIDs( Symmetry, idSymmetry )

idInflatedBubble = []
InflatedBubble = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idInflatedBubble.append( returnIDofShape( [ 0,     h_s, 0], "EDGE" ) )
InflatedBubble_union = geompy.UnionIDs( InflatedBubble, idInflatedBubble )


idAmbient = []
Ambient = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idAmbient.append( returnIDofShape( [    -Height_needle, Radius_needle+h_s, 0] , "EDGE" ) )
idAmbient.append( returnIDofShape( [    -Height_needle, Radius_needle+dR_ref+h_s, 0] , "EDGE" ) )
idAmbient.append( returnIDofShape( [    -Height_needle,   Radius_tank-h_s, 0] , "EDGE" ) )
# idAmbient.append( returnIDofShape( [ Height_tank,     h_s, 0] , "EDGE" ) )
Ambient_union = geompy.UnionIDs( Ambient, idAmbient )

idInternalEquidistributionZ = []
InternalEquidistributionZ = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idInternalEquidistributionZ.append( returnIDofShape( [ h_s, Radius_needle, 0] , "EDGE" ) )
InternalEquidistributionZ_union = geompy.UnionIDs( InternalEquidistributionZ, idInternalEquidistributionZ)

idInternalEquidistributionR = []
InternalEquidistributionR = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idInternalEquidistributionR.append( returnIDofShape( [ 0, Radius_needle+h_s, 0] , "EDGE" ) )
InternalEquidistributionR_union = geompy.UnionIDs( InternalEquidistributionR, idInternalEquidistributionR)


idverticals1 = []
verticals1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idverticals1.append( idAmbient[0] )
idverticals1.append( idInternalEquidistributionR[0] )
idverticals1.append(  returnIDofShape( [dR_ref,        Radius_needle+h_s, 0], "EDGE" ) )
verticals1_union = geompy.UnionIDs( verticals1, idverticals1 )



idverticals2 = []
verticals2 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idverticals2.append( idInflatedBubble[0] )
idverticals2.append( returnIDofShape( [ dR_ref,        h_s, 0 ] , "EDGE") )
verticals2_union = geompy.UnionIDs( verticals2, idverticals2 )

idverticals3 = []
verticals3 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idverticals3.append( returnIDofShape( [ 0, Radius_needle-Needle_thickness+h_s, 0 ] , "EDGE") )
idverticals3.append( returnIDofShape( [ dR_ref, Radius_needle-Needle_thickness+h_s, 0 ] , "EDGE") )
verticals2_union = geompy.UnionIDs( verticals3, idverticals3 )

geompy.addToStudy(verticals3,'b')


idhorizontals1 = []
horizontals1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idhorizontals1.append( idWall[0] )
idhorizontals1.append( returnIDofShape( [-h_s,        Radius_needle + dR_ref, 0], "EDGE" ) )
horizontals1_union = geompy.UnionIDs( horizontals1, idhorizontals1 )


idhorizontals2 = []
horizontals2 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idhorizontals2.append( returnIDofShape( [ h_s, Radius_needle + dR_ref, 0] , "EDGE" ) )
idhorizontals2.append( returnIDofShape( [ h_s,        Radius_needle, 0], "EDGE" ) )
idhorizontals2.append( returnIDofShape( [ h_s, Radius_needle-Needle_thickness, 0], "EDGE" ) )
idhorizontals2.append( idSymmetry[0] )
horizontals2_union = geompy.UnionIDs( horizontals2, idhorizontals2 )

geompy.addToStudy(horizontals2,'a')
geompy.addToStudy(verticals2,'a')

idquad1 = []
quad1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idquad1.append( returnIDofShape( [ -h_s, Radius_needle+h_s, 0] , "FACE" ) )
quad1_union = geompy.UnionIDs( quad1, idquad1 )

idquad2 = []
quad2 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idquad2.append( returnIDofShape( [ +h_s, Radius_needle+h_s, 0] , "FACE" ) )
quad2_union = geompy.UnionIDs( quad2, idquad2 )


idquad3 = []
quad3 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idquad3.append( returnIDofShape( [ +h_s, +h_s, 0] , "FACE" ) )
quad3_union = geompy.UnionIDs( quad3, idquad3 )


idquad4 = []
quad4 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idquad4.append( returnIDofShape( [ +h_s, Radius_needle-Needle_thickness+h_s, 0] , "FACE" ) )
quad3_union = geompy.UnionIDs( quad4, idquad4 )
geompy.addToStudy(quad4,'c')

idsquare = []
square = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idsquare.append( returnIDofShape( [    -Height_needle + h_s, Radius_needle+dR_ref+h_s, 0] , "FACE" ) )
square_union = geompy.UnionIDs( square, idsquare )



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


Mesh_1.Segment(geom=verticals2).NumberOfSegments(NumSegmentsOnInflatedBubble)
# Mesh_1.Segment(geom=verticals1).StartEndLength (startLengthOnSurroundings,endLengthOnSurroundings,[])
Mesh_1.Segment(geom=verticals1).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings,[])
Mesh_1.Segment(geom=verticals3).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings,[])
Mesh_1.Segment(geom=horizontals1).StartEndLength (startLengthOnSurroundings,endLengthOnSurroundings, idhorizontals1)
Mesh_1.Segment(geom=horizontals2).NumberOfSegments(NumSegmentsOnSymmetry)

SM4 = Mesh_1.Quadrangle(geom = quad1)
SM5 = Mesh_1.Quadrangle(geom = quad2)
SM7 = Mesh_1.Quadrangle(geom = quad3)
SM7 = Mesh_1.Quadrangle(geom = quad4)

NETGEN_1D_2D_global   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
SM6   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = square)

MeshParameters(NETGEN_1D_2D_global  ,Main_maxSize_element, Main_minSize_element, 0.2)
MeshParameters(SM6  ,maxSquareElementLength, minSquareElementLength, 0.3)


GSM4 = SM4.GetSubMesh()
GSM5 = SM5.GetSubMesh()
GSM6 = SM6.GetSubMesh()
GSM7 = SM7.GetSubMesh()


isDone = Mesh_1.SetMeshOrder( [   [ GSM4,
                                    GSM5,
                                    GSM7,
                                    GSM6] ])

Wall_1           			=  Mesh_1.GroupOnGeom( Wall 	,'Wall' 	 ,SMESH.EDGE )
Symmetry_1       			=  Mesh_1.GroupOnGeom( Symmetry ,'Symmetry'	 ,SMESH.EDGE )
InflatedBubble_1 			=  Mesh_1.GroupOnGeom( InflatedBubble 	,'InflatedBubble' 	 ,SMESH.EDGE )
Ambient_1        			=  Mesh_1.GroupOnGeom( Ambient 	,'Ambient' 	 ,SMESH.EDGE )
InternalEquidistributionR_1 =  Mesh_1.GroupOnGeom( InternalEquidistributionR 	,'InternalEq_R' 	 ,SMESH.EDGE )
InternalEquidistributionZ_1 =  Mesh_1.GroupOnGeom( InternalEquidistributionZ 	,'InternalEq_Z' 	 ,SMESH.EDGE )


isDone = Mesh_1.Compute()


isDone = Mesh_1.SplitQuadObject( Mesh_1, 1 )

# isDone = Mesh_1.QuadTo4Tri( )
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
