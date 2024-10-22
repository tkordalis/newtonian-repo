#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.5.0 with dump python functionality
###

import sys
import salome 
import read_datfile


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

from Geometry_Mesh_Parameters import Radius_tank, Height_tank, RSphere1, dR_ref1, dR_ref2,  dR_ref3, dR_ref4, R_refinement1_Sphere1, R_refinement2_Sphere1, R_refinement3_Sphere1, R_refinement4_Sphere1, \
ellipse_position, ellipse_Minor_Radius, ellipse_Major_Radius, outer_ellipse_Major_Radius, outer_ellipse_Minor_Radius, h_s, \
Main_maxSize_element, Main_minSize_element, Element_size_on_Sphere, Netgen_Params, NumSegmentsOnSphere, Element_size_on_Ambient, NodeDensityFunction_Sym, NumSegmentsOnAmbient, \
dZ_refAmb, Element_size_on_Ambient_cb

def min_distance_index_value(p_o, p_rest):
    # creates a list of the distances of all point to our reference point
    # and returns the index of the min value
    min_dist_index_value = []
    dist = []
    for item in p_rest:
        dist.append(hypot(item[0]-p_o[0], item[1]-p_o[1]))

    dist_temp = sorted(dist)
    min_dist = dist_temp[0]
    for i in range(len(dist)):
        if abs(min_dist - dist[i])<0.000001:
            min_dist_index = i
            break
    min_dist_index_value.append(min_dist_index)
    min_dist_index_value.append(min_dist)

    return min_dist_index_value


def sortCoordinatesOfBoundary(CoordOfBoundary, VariableToSort):
    # coordinates of point with R = 0
    firstval = 0
    if VariableToSort == 1:
        maxZ_point_coords = max(CoordOfBoundary, key=lambda tup: tup[0])
    elif VariableToSort == 2:
        maxZ_point_coords = max(CoordOfBoundary, key=lambda tup: tup[1])


    firstval = maxZ_point_coords
    

    
    nextval = firstval
    returnList = []
    CoordOfBoundary_temp  = []
    CoordOfBoundary_temp  = CoordOfBoundary
    CoordOfBoundarySorted = []

    ArcLength = 0.
    for iterate in range(len(CoordOfBoundary)-1):
        CoordOfBoundary_temp.remove(nextval)
        CoordOfBoundarySorted.append(nextval)
        index_next_point = min_distance_index_value(nextval, CoordOfBoundary_temp)[0]
        dS = min_distance_index_value(nextval, CoordOfBoundary_temp)[1]
        
        ArcLength = ArcLength + dS
        nextval = CoordOfBoundary_temp[index_next_point]

    CoordOfBoundarySorted.append(CoordOfBoundary_temp[0])

    returnList.append(CoordOfBoundarySorted)
    returnList.append(ArcLength)
    # print(ArcLength)
    return returnList


def CurveCoordinatesScaling(CoordsOfCurve_input, DistanceFactor):

    ScaledCurveCoordinates = []

    CoordsOfCurve = CoordsOfCurve_input.copy()
    # CoordsOfCurve.reverse()
    point_start_coords = CoordsOfCurve[0]
    # i put minus because i start the points from the negative
    ScaledCurveCoordinates.append( [ point_start_coords[0] - DistanceFactor, point_start_coords[1] ] )

    for i in range( len(CoordsOfCurve) - 1 ):

        midpoint_coords     = []
        tangent_vector      = []
        normal_vector       = []
        unit_normal_vector  = []
        ScaleFactor         = []
        scaled_point_coords = []


        point1coords = CoordsOfCurve[i]
        point2coords = CoordsOfCurve[i+1]

        x_coord_point_mid = 0.5* ( point2coords[0] + point1coords[0] )
        y_coord_point_mid = 0.5* ( point2coords[1] + point1coords[1] )

        midpoint_coords.append( x_coord_point_mid )
        midpoint_coords.append( y_coord_point_mid )

        # tangent_vector.append( point2coords[0] - point1coords[0] )
        # tangent_vector.append( point2coords[1] - point1coords[1] )
        tangent_vector.append( point1coords[0] - point2coords[0] )
        tangent_vector.append( point1coords[1] - point2coords[1] )

        normal_vector.append(  tangent_vector[1] )
        normal_vector.append( -tangent_vector[0] )
        unit_normal_vector.append( normal_vector[0]/ sqrt( normal_vector[0]**2 + normal_vector[1]**2 ) )
        unit_normal_vector.append( normal_vector[1]/ sqrt( normal_vector[0]**2 + normal_vector[1]**2 ) )
        # print(unit_normal_vector)
        ScaleFactor.append( DistanceFactor*unit_normal_vector[0] )
        ScaleFactor.append( DistanceFactor*unit_normal_vector[1] )

        scaled_point_coords.append( ScaleFactor[0]+midpoint_coords[0] )
        scaled_point_coords.append( ScaleFactor[1]+midpoint_coords[1] )

        ScaledCurveCoordinates.append(scaled_point_coords)



    point_end_coords = CoordsOfCurve[-1]

    ScaledCurveCoordinates.append( [ point_end_coords[0]+ DistanceFactor, point_end_coords[1] ] )

    return ScaledCurveCoordinates

geompy = geomBuilder.New()

segment_length_from_NumberOfSegments = 3.14159265359/NumSegmentsOnSphere
cb1NumSegments = 6*int(dR_ref1/segment_length_from_NumberOfSegments)


# Base Vectors
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(7, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 7, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 7)


Bubble1PointCoordinates = []
for p in read_datfile.readBoundaryNodes("Bubble1.dat"):
    Bubble1PointCoordinates.append(p)

Bubble2PointCoordinates = []
for p in read_datfile.readBoundaryNodes("Bubble2.dat"):
    Bubble2PointCoordinates.append(p)


AmbientPointCoordinates = []
for p in read_datfile.readBoundaryNodes("Ambient.dat"):
    AmbientPointCoordinates.append(p)


result_B1 = sortCoordinatesOfBoundary(Bubble1PointCoordinates, 1)
Bubble1PointCoordinatesSorted = result_B1[0]

result_B2 = sortCoordinatesOfBoundary(Bubble2PointCoordinates, 1)
Bubble2PointCoordinatesSorted = result_B2[0]

result_Amb = sortCoordinatesOfBoundary(AmbientPointCoordinates, 2)
AmbientPointCoordinatesSorted = result_Amb[0]


Bubble1Points = [geompy.MakeVertex(*p, 0) for p in Bubble1PointCoordinatesSorted]
Bubble2Points = [geompy.MakeVertex(*p, 0) for p in Bubble2PointCoordinatesSorted]
i=-1
for item in Bubble1Points:
    i+=1
    geompy.addToStudy(item,str(i))


B1_leftPoint  = Bubble1Points[-1]
B1_rightPoint = Bubble1Points[0]

B2_leftPoint  = Bubble2Points[-1]
B2_rightPoint = Bubble2Points[0]

distance = (Bubble1PointCoordinatesSorted[0])[0]-(Bubble2PointCoordinatesSorted[-1])[0]

AmbientPoints = [geompy.MakeVertex(*p, 0) for p in AmbientPointCoordinatesSorted]


geompy.addToStudy(B1_leftPoint,'B1_leftPoint')

Curve_1 = geompy.MakeInterpol(Bubble1Points, False)
B1_line = geompy.MakeLineTwoPnt(B1_leftPoint, B1_rightPoint)
B1_Face = geompy.MakeFaceWires ([Curve_1, B1_line], 1)

geompy.addToStudy(B1_Face,'B1_Face')

Curve_2 = geompy.MakeInterpol(Bubble2Points, False)
B2_line = geompy.MakeLineTwoPnt(B2_leftPoint, B2_rightPoint)
B2_Face = geompy.MakeFaceWires ([Curve_2, B2_line], 1)


bottomRightVertexTank = geompy.MakeVertex( *AmbientPointCoordinatesSorted[-1] , 0 )
topRightVertexTank    = geompy.MakeVertex( *AmbientPointCoordinatesSorted[0], 0 )

bottomLeftVertexTank  = geompy.MakeVertex( -0.5*Height_tank, 0, 0 )
topLeftVertexTank  = geompy.MakeVertex( -0.5*Height_tank, Radius_tank, 0 )

bottomLineTank = geompy.MakeLineTwoPnt( bottomLeftVertexTank , bottomRightVertexTank )
# rightLineTank  = geompy.MakeLineTwoPnt( bottomRightVertexTank, topRightVertexTank )
rightLineTank  = geompy.MakeInterpol(AmbientPoints, False)

topLineTank   = geompy.MakeLineTwoPnt( topRightVertexTank, topLeftVertexTank )
leftLineTank  = geompy.MakeLineTwoPnt( topLeftVertexTank , bottomLeftVertexTank )

Domain =  geompy.MakeFaceWires([ geompy.MakeWire([bottomLineTank, rightLineTank, topLineTank, leftLineTank], 1e-08) ] , 1)

geompy.addToStudy(Domain,'Domain')

# halfCut = geompy.MakeTranslation( geompy.MakeFaceHW(Height_tank, 2*Radius_tank, 1), 0, -Radius_tank, 0 )
Domain_cut = geompy.MakeCutList( Domain, [ B1_Face, B2_Face ], True )

# geompy.addToStudy(Domain,'Domain')
# geompy.addToStudy(Domain_cut,'Domain_cut')

Scaled_Curve_Bubble1_Ref1_coords = []
Scaled_Curve_Bubble1_Ref1_points = []
Scaled_Curve_Bubble1_Ref1_coords = CurveCoordinatesScaling(Bubble1PointCoordinatesSorted, -dR_ref1)
Scaled_Curve_Bubble1_Ref1_points = [geompy.MakeVertex(*p, 0) for p in Scaled_Curve_Bubble1_Ref1_coords]
Scaled_Curve_Bubble1_Ref1 = geompy.MakeInterpol(Scaled_Curve_Bubble1_Ref1_points, False)

Scaled_Curve_Bubble1_Ref2_coords = []
Scaled_Curve_Bubble1_Ref2_points = []
Scaled_Curve_Bubble1_Ref2_coords = CurveCoordinatesScaling(Bubble1PointCoordinatesSorted, -dR_ref2)
Scaled_Curve_Bubble1_Ref2_points = [geompy.MakeVertex(*p, 0) for p in Scaled_Curve_Bubble1_Ref2_coords]
Scaled_Curve_Bubble1_Ref2 = geompy.MakeInterpol(Scaled_Curve_Bubble1_Ref2_points, False)


Scaled_Curve_Bubble2_Ref1_coords = []
Scaled_Curve_Bubble2_Ref1_points = []
Scaled_Curve_Bubble2_Ref1_coords = CurveCoordinatesScaling(Bubble2PointCoordinatesSorted, -dR_ref1)
Scaled_Curve_Bubble2_Ref1_points = [geompy.MakeVertex(*p, 0) for p in Scaled_Curve_Bubble2_Ref1_coords]
Scaled_Curve_Bubble2_Ref1 = geompy.MakeInterpol(Scaled_Curve_Bubble2_Ref1_points, False)

Scaled_Curve_Bubble2_Ref2_coords = []
Scaled_Curve_Bubble2_Ref2_points = []
Scaled_Curve_Bubble2_Ref2_coords = CurveCoordinatesScaling(Bubble2PointCoordinatesSorted, -dR_ref2)
Scaled_Curve_Bubble2_Ref2_points = [geompy.MakeVertex(*p, 0) for p in Scaled_Curve_Bubble2_Ref2_coords]
Scaled_Curve_Bubble2_Ref2 = geompy.MakeInterpol(Scaled_Curve_Bubble2_Ref2_points, False)

R_refinement3_Sphere1 = 0.5*distance + dR_ref3
R_refinement4_Sphere1 = 0.5*distance + dR_ref4
midpoint = (Bubble2PointCoordinatesSorted[-1])[0] + 0.5*distance
Disk_refinement3_Sphere1 = geompy.MakeTranslation( geompy.MakeDiskR(R_refinement3_Sphere1,1), midpoint, 0, 0 )
Disk_refinement4_Sphere1 = geompy.MakeTranslation( geompy.MakeDiskR(R_refinement4_Sphere1,1), midpoint, 0, 0 )
[Wire_3] = geompy.ExtractShapes(Disk_refinement3_Sphere1, geompy.ShapeType["WIRE"], True)
[Wire_4] = geompy.ExtractShapes(Disk_refinement4_Sphere1, geompy.ShapeType["WIRE"], True)


# Scaled_Curve_Bubble1_Ref3_coords = []
# Scaled_Curve_Bubble1_Ref3_points = []
# Scaled_Curve_Bubble1_Ref3_coords = CurveCoordinatesScaling(Bubble1PointCoordinatesSorted, -dR_ref3)
# Scaled_Curve_Bubble1_Ref3_points = [geompy.MakeVertex(*p, 0) for p in Scaled_Curve_Bubble1_Ref3_coords]
# Scaled_Curve_Bubble1_Ref3 = geompy.MakeInterpol(Scaled_Curve_Bubble1_Ref3_points, False)

# Scaled_Curve_Bubble1_Ref4_coords = []
# Scaled_Curve_Bubble1_Ref4_points = []
# Scaled_Curve_Bubble1_Ref4_coords = CurveCoordinatesScaling(Bubble1PointCoordinatesSorted, -dR_ref4)
# Scaled_Curve_Bubble1_Ref4_points = [geompy.MakeVertex(*p, 0) for p in Scaled_Curve_Bubble1_Ref4_coords]
# Scaled_Curve_Bubble1_Ref4 = geompy.MakeInterpol(Scaled_Curve_Bubble1_Ref4_points, False)

geompy.addToStudy(Scaled_Curve_Bubble1_Ref1,'Scaled_Curve_Bubble1_Ref1')
geompy.addToStudy(Scaled_Curve_Bubble1_Ref2,'Scaled_Curve_Bubble1_Ref2')
geompy.addToStudy(Scaled_Curve_Bubble2_Ref1,'Scaled_Curve_Bubble1_Ref1')
geompy.addToStudy(Scaled_Curve_Bubble2_Ref2,'Scaled_Curve_Bubble1_Ref2')


Centroid        = geompy.MakeVertexOnCurve(B1_line, 0.5, True)
Centroid_coords = geompy.PointCoordinates(Centroid)



partition_line_list = [ Scaled_Curve_Bubble1_Ref1, Scaled_Curve_Bubble1_Ref2, Scaled_Curve_Bubble2_Ref1, Scaled_Curve_Bubble2_Ref2, Wire_3, Wire_4 ]

partition_tool = geompy.MakeFuseList( partition_line_list, True, True)

Partition_1 = geompy.MakePartition([Domain_cut], [partition_tool], [], [], geompy.ShapeType["FACE"], 0, [], 0)

geompy.addToStudy(Partition_1,'Partition_1')


def returnPointsFromRotatedAxis( theta_degrees, x_tilt, x_o ):
    x_global        = []
    theta_radians   = 3.1415926*(theta_degrees/180.0)

    x_global.append(x_o[0] + cos(theta_radians)*x_tilt[0] + (-sin(theta_radians))*x_tilt[1])
    x_global.append(x_o[1] + sin(theta_radians)*x_tilt[0] +   cos(theta_radians) *x_tilt[1])
    x_global.append(0.0)

    return x_global

def returnIDofShape( theta_degrees, x_tilt, x_o, TypeofShape ):
    idShape              =  0.
    Shape_Point_Position =  returnPointsFromRotatedAxis( theta_degrees, x_tilt, x_o )
    Shape_Point          =  geompy.MakeVertex(*Shape_Point_Position)
    sLineShape           =  geompy.GetShapesNearPoint(Partition_1, Shape_Point, geompy.ShapeType[ TypeofShape ])
    idShape              =  geompy.GetSubShapeID(Partition_1,   sLineShape  )
    # a1 = Shape_Point
    # geompy.addToStudy(a1, "a1")

    return idShape

def getgroupSymmetry( mainGroup, refZone, leftOrRight, gid ):
    groupSymmetryDict = {
    "mainGroup"   : mainGroup,
    "refZone"     : refZone,
    "leftOrRight" : leftOrRight,
    "id"          : gid
    }
    return groupSymmetryDict

def getgroupFace( mainGroup, refZone, gid, obj, union, mesh_obj):
    groupSphereDict = {
    "mainGroup"   : mainGroup,
    "refZone"       : refZone,
    "id"                : gid,
    "obj"               : obj,
    "union"         : union,
    "mesh_obj"      : mesh_obj
    }
    return groupSphereDict


idtankWall = []
tankWall       = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idtankWall.append(returnIDofShape(0, [ 0, Radius_tank ], [0, 0], "EDGE"))
idtankWall.append( returnIDofShape(0, [ -0.5*Height_tank, +h_s ], [0, 0], "EDGE") ) 
idtankWall.append( returnIDofShape( 0, [ 0.5*Height_tank-h_s,Radius_tank ], [0, 0], "EDGE" ) )
tankWall_union = geompy.UnionIDs(tankWall , idtankWall  )


idrightPlane = [] 
rightPlane   = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idrightPlane.append( returnIDofShape(0, [ (AmbientPointCoordinates[0])[0], +h_s ], [0, 0], "EDGE") )
rightPlane_union = geompy.UnionIDs( rightPlane , idrightPlane )

idSphere1     = []
idSphere1.append(returnIDofShape( 0, [(geompy.PointCoordinates(Bubble1Points[5]))[0], (geompy.PointCoordinates(Bubble1Points[5]))[1]], [0,0], "EDGE" ))
Sphere1       = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Sphere1_union = geompy.UnionIDs( Sphere1 , idSphere1 )

idSphere2     = []
idSphere2.append(returnIDofShape( 0, [(geompy.PointCoordinates(Bubble2Points[5]))[0], (geompy.PointCoordinates(Bubble2Points[5]))[1]], [0,0], "EDGE" ))
Sphere2       = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Sphere2_union = geompy.UnionIDs( Sphere2 , idSphere2 )


idSymmetry = []

Symmetry_groups = []


Symmetry_groups.append( getgroupSymmetry( 'middle', None     , None , None) )

Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb1'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb1'     , 'right', None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb2'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'cb2'     , 'right', None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'ellipse1', 'right' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'ellipse2', 'right' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_1_Side', 'out'     , 'right' , None) )

Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'cb1'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'cb1'     , 'right', None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'cb2'     , 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'cb2'     , 'right', None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'ellipse1', 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'ellipse2', 'left' , None) )
Symmetry_groups.append( getgroupSymmetry( 'Sphere_2_Side', 'out'     , 'left' , None) )



origin_of_axes = []
x_tilt         = []
x_tilt_dumy    = []
theta_degrees  = []

point_list = []
x_tilt = [0,0]

Symmetry_groups[0]["id"] = returnIDofShape( 0, [midpoint,0] , [0,0], "EDGE" )
idSymmetry.append( Symmetry_groups[0]["id"] )
for i in range(1,len(Symmetry_groups)):

    origin_of_axes = [0,0]
    theta_degrees = 0
    if Symmetry_groups[i]["mainGroup"] == "Sphere_1_Side":

        if Symmetry_groups[i]["refZone"] == "cb1" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  Bubble1PointCoordinatesSorted[-1]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]

        elif Symmetry_groups[i]["refZone"] == "cb1" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  Bubble1PointCoordinatesSorted[0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "cb2" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  Scaled_Curve_Bubble1_Ref1_coords[-1]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "cb2" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  Scaled_Curve_Bubble1_Ref1_coords[0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "ellipse1" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  Scaled_Curve_Bubble1_Ref2_coords[0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "ellipse2" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  [midpoint + R_refinement3_Sphere1, 0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "out" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  [midpoint + R_refinement4_Sphere1, 0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]

    elif Symmetry_groups[i]["mainGroup"] == "Sphere_2_Side":

        if Symmetry_groups[i]["refZone"] == "cb1" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  Bubble2PointCoordinatesSorted[-1]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]

        elif Symmetry_groups[i]["refZone"] == "cb1" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  Bubble2PointCoordinatesSorted[0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "cb2" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  Scaled_Curve_Bubble2_Ref1_coords[-1]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "cb2" and Symmetry_groups[i]["leftOrRight"] == "right":
            x_tilt_dumy =  Scaled_Curve_Bubble2_Ref1_coords[0]
            x_tilt[0] =  x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "ellipse1" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  Scaled_Curve_Bubble2_Ref2_coords[-1]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "ellipse2" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  [midpoint - R_refinement3_Sphere1, 0]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]


        elif Symmetry_groups[i]["refZone"] == "out" and Symmetry_groups[i]["leftOrRight"] == "left":
            x_tilt_dumy =  [midpoint - R_refinement4_Sphere1, 0]
            x_tilt[0] =  x_tilt_dumy[0]-h_s
            x_tilt[1] = x_tilt_dumy[1]
            
    # if Symmetry_groups[i]["leftOrRight"] == "left":
    #     theta_degrees = 180
    # elif Symmetry_groups[i]["leftOrRight"] == "right":
    #     theta_degrees = 0

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

geompy.addToStudyInFather(Partition_1, Symmetry, "Symmetry")




idSphere1Ref1     = [returnIDofShape( 0, Scaled_Curve_Bubble1_Ref1_coords[5], [0,0], "EDGE" )]
Sphere1Ref1         = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Sphere1_unionRef1 = geompy.UnionIDs( Sphere1Ref1 , idSphere1Ref1 )

idSphere2Ref1     = [returnIDofShape( 0, Scaled_Curve_Bubble2_Ref1_coords[5], [0,0], "EDGE" )]
Sphere2Ref1         = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
Sphere2_unionRef1 = geompy.UnionIDs( Sphere2Ref1 , idSphere2Ref1 )

geompy.addToStudyInFather(Partition_1,Sphere1,"Sphere1")
geompy.addToStudyInFather(Partition_1,Sphere1Ref1,"Sphere1Ref1")

geompy.addToStudyInFather(Partition_1,Sphere2,"Sphere2")
geompy.addToStudyInFather(Partition_1,Sphere2Ref1,"Sphere1Ref2")


Groups_faces = []

Groups_faces.append( getgroupFace("Sphere1", "Ref_1", None, None, None, None) )
Groups_faces.append( getgroupFace("Sphere2", "Ref_1", None, None, None, None) )
Groups_faces.append( getgroupFace("Sphere1", "Ref_2", None, None, None, None) )
Groups_faces.append( getgroupFace("Sphere2", "Ref_2", None, None, None, None) )
Groups_faces.append( getgroupFace("ellipse", "1"    , None, None, None, None) )
Groups_faces.append( getgroupFace("ellipse", "2"    , None, None, None, None) )

for i in range(len(Groups_faces)):

    if Groups_faces[i]["mainGroup"] == "Sphere1":
        origin_of_axes = [0,0]

        if Groups_faces[i]["refZone"] == "Ref_1":
            x_tilt_dumy = Bubble1PointCoordinatesSorted[5]
            x_tilt[0] = x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]+h_s
        elif Groups_faces[i]["refZone"] == "Ref_2":
            x_tilt_dumy = Scaled_Curve_Bubble1_Ref1_coords[5]
            x_tilt[0] = x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]+h_s
    elif Groups_faces[i]["mainGroup"] == "Sphere2":
        if Groups_faces[i]["refZone"] == "Ref_1":
            x_tilt_dumy = Bubble2PointCoordinatesSorted[5]
            x_tilt[0] = x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]+h_s
        elif Groups_faces[i]["refZone"] == "Ref_2":
            x_tilt_dumy = Scaled_Curve_Bubble2_Ref1_coords[5]
            x_tilt[0] = x_tilt_dumy[0]+h_s
            x_tilt[1] = x_tilt_dumy[1]+h_s
    elif Groups_faces[i]["mainGroup"] == "ellipse":
        if Groups_faces[i]["refZone"] == "1":
            x_tilt_dumy = [midpoint, R_refinement3_Sphere1-h_s]
            x_tilt[0] = x_tilt_dumy[0]
            x_tilt[1] = x_tilt_dumy[1]
        elif Groups_faces[i]["refZone"] == "2":
            x_tilt_dumy = [midpoint, R_refinement4_Sphere1-h_s]
            x_tilt[0] = x_tilt_dumy[0]
            x_tilt[1] = x_tilt_dumy[1]


    theta_degrees = 0

    Groups_faces[i]["id"] = returnIDofShape( theta_degrees, x_tilt, origin_of_axes, "FACE")

i=-1
for item in Groups_faces:
    i = i+1
    item["obj"] = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
    item["union"] = geompy.UnionIDs(item["obj"] , [ item["id" ]] )
    # print(item["mainGroup"], item["refZone"])
    geompy.addToStudyInFather(Partition_1, item["obj"], str(i))
    # print(Groups_faces[i]["id"])


# # # --------------------------- End of Geometry --------------------------- #

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
Mesh_1.Segment(geom=Sphere1).NumberOfSegments(NumSegmentsOnSphere)
Mesh_1.Segment(geom=Sphere2).NumberOfSegments(NumSegmentsOnSphere)

Mesh_1.Segment(geom=Sphere1Ref1).NumberOfSegments(NumSegmentsOnSphere)
Mesh_1.Segment(geom=Sphere2Ref1).NumberOfSegments(NumSegmentsOnSphere)
Mesh_1.Segment(geom=horizontalsSphere1).NumberOfSegments(cb1NumSegments)
NETGEN_1D_2D   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
MeshParameters(NETGEN_1D_2D  ,Main_maxSize_element, Main_minSize_element, 0.1)

Params = []

Groups_faces[0]["mesh_obj"] = Mesh_1.Quadrangle(geom = Groups_faces[0]["obj"])
Groups_faces[1]["mesh_obj"] = Mesh_1.Quadrangle(geom = Groups_faces[1]["obj"]) 
j=0
for i in range(3,len(Groups_faces)):
    # print(str(i))
    Params = Netgen_Params[j]
    Groups_faces[i]["mesh_obj"] = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = Groups_faces[i]["obj"])
    MeshParameters(Groups_faces[i]["mesh_obj"], Params[0], Params[1], Params[2])
    j+=1

Params = Netgen_Params[0]
Groups_faces[2]["mesh_obj"] = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = Groups_faces[2]["obj"])
MeshParameters(Groups_faces[2]["mesh_obj"], Params[0], Params[1], Params[2])

SymmOut_params = Netgen_Params[-1]
Mesh_1.Segment(geom=SymmetryB1out).StartEndLength ( Main_maxSize_element, SymmOut_params[0] )
Mesh_1.Segment(geom=SymmetryB2out).StartEndLength ( SymmOut_params[0], Main_maxSize_element )

tankWall_1      =  Mesh_1.GroupOnGeom( tankWall     ,'tankWall'   ,SMESH.EDGE )
Symmetry_1      =  Mesh_1.GroupOnGeom( Symmetry     ,'Symmetry'   ,SMESH.EDGE )
Sphere1_1       =  Mesh_1.GroupOnGeom( Sphere1      ,'Bubble1'    ,SMESH.EDGE )
Sphere2_1       =  Mesh_1.GroupOnGeom( Sphere2      ,'Bubble2'    ,SMESH.EDGE )
rightPlane_1    =  Mesh_1.GroupOnGeom( rightPlane ,'Ambient'        ,SMESH.EDGE )


Priority_list = []

for item in Groups_faces:
    SubMesh = item["mesh_obj"]
    # print(SubMesh)
    Priority_list.append(SubMesh.GetSubMesh())

# print(Priority_list)
isDone = Mesh_1.SetMeshOrder( [Priority_list] )

isDone = Mesh_1.Compute()

# isDone = Mesh_1.QuadTo4Tri( )
isDone = Mesh_1.SplitQuadObject( Mesh_1, 0 )



try:
  Mesh_1.ExportUNV( r'./Bounded.unv' )
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
