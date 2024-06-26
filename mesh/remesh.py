#####
#
#
###
### This file is generated automatically by SALOME v9.5.0 with dump python functionality
###

import sys
import salome
import read_datfile

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/tkordalis')

###
### GEOM component
###
import GEOM
from salome.geom import geomBuilder
from math import *
import SALOMEDS
from Geometry_Mesh_Parameters import *



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



def sortCoordinatesOfBoundary(CoordOfBoundary):
    # coordinates of point with R = 0
    firstval = 0
    minR_point_coords = min(CoordOfBoundary, key=lambda tup: tup[1])
    firstval = minR_point_coords
    # for item in CoordOfBoundary:
    #     if item[1] < 0.000001:
    #         firstval = item
    #         break

    
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
    return returnList


def CurveCoordinatesScaling(CoordsOfCurve, DistanceFactor):

    ScaledCurveCoordinates = []

    point_start_coords = CoordsOfCurve[0]


    # point_start_coords[0] = point_start_coords[0] + DistanceFactor
    # point_start_coords[1] = point_start_coords[1]


    ScaledCurveCoordinates.append( [ point_start_coords[0] + DistanceFactor, point_start_coords[1] ] )

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

        tangent_vector.append( point2coords[0] - point1coords[0] )
        tangent_vector.append( point2coords[1] - point1coords[1] )

        normal_vector.append(  tangent_vector[1] )
        normal_vector.append( -tangent_vector[0] )

        unit_normal_vector.append( normal_vector[0]/ sqrt( normal_vector[0]**2 + normal_vector[1]**2 ) )
        unit_normal_vector.append( normal_vector[1]/ sqrt( normal_vector[0]**2 + normal_vector[1]**2 ) )
        ScaleFactor.append( DistanceFactor*unit_normal_vector[0] )
        ScaleFactor.append( DistanceFactor*unit_normal_vector[1] )

        scaled_point_coords.append( ScaleFactor[0]+midpoint_coords[0] )
        scaled_point_coords.append( ScaleFactor[1]+midpoint_coords[1] )

        ScaledCurveCoordinates.append(scaled_point_coords)



    point_end_coords = CoordsOfCurve[-1]

    ScaledCurveCoordinates.append( [ point_end_coords[0]- DistanceFactor, point_end_coords[1] ] )

    return ScaledCurveCoordinates









# there will be one big if 
# if the radial coordinate of the highest point of the free surface is smaller than the  radius of the needle then 
# the mesh is the only i have already constructed
# BUT if the named radial coordinate is LARGER than the radius of the needle then i will construct only one structured area around the 
# free surface with partition tool ONLY the !!-Total_Scaled_Curve_Bubble1_Ref1-!! 


# global stuff

Bubble1PointCoordinates = []
for p in read_datfile.readBoundaryNodes("Bubble1.dat"):
    Bubble1PointCoordinates.append(p)
    # print(p)

Bubble1PointCoordinatesSorted = sortCoordinatesOfBoundary(Bubble1PointCoordinates)[0]
B1_ArcLength = sortCoordinatesOfBoundary(Bubble1PointCoordinates)[1]



# for item in Bubble1PointCoordinatesSorted:
#     print(item)

AmbientPointCoordinates = []
for p in read_datfile.readBoundaryNodes("Ambient.dat"):
    AmbientPointCoordinates.append(p)
    

dz_ambient = abs( ( AmbientPointCoordinates[0] )[0] + Height_needle )
print(dz_ambient)
geompy = geomBuilder.New()

Height_tank_modified   = Height_tank + dz_ambient
Height_needle_modified = Height_needle + dz_ambient

test1 = geompy.MakeFaceHW(Height_tank, Radius_tank, 1)
test2 = geompy.MakeFaceHW(Height_tank_modified, Radius_tank, 1)
geompy.addToStudy(test1,'test1')
geompy.addToStudy(test2,'test2')

# Domain_cut_uninflated = geompy.MakeCutList( geompy.MakeTranslation( geompy.MakeFaceHW(Height_tank, Radius_tank, 1), 0.5*Height_tank-Height_needle, 0.5*Radius_tank, 0 ), \
# 								[ geompy.MakeTranslation( geompy.MakeFaceHW(Height_needle, Radius_needle, 1), -0.5*Height_needle, 0.5*Radius_needle, 0 ) ], True )

Domain_cut_uninflated = geompy.MakeCutList( geompy.MakeTranslation( geompy.MakeFaceHW(Height_tank_modified, Radius_tank, 1), 0.5*Height_tank-Height_needle-0.5*dz_ambient, 0.5*Radius_tank, 0 ), \
                                [ geompy.MakeTranslation( geompy.MakeFaceHW(Height_needle_modified, Radius_needle, 1), -0.5*Height_needle_modified, 0.5*Radius_needle, 0 ) ], True )
geompy.addToStudy(Domain_cut_uninflated,'Domain_cut_uninflated')

# geompy.addToStudy(Domain_cut_uninflated2, 'Domain_cut_uninflated2')

Bubble1Points = [geompy.MakeVertex(*p, 0) for p in Bubble1PointCoordinatesSorted]

Curve_1 = geompy.MakeInterpol(Bubble1Points, False)
geompy.addToStudy(Curve_1,'Curve_1')

horizontal_line = geompy.MakeLineTwoPnt( geompy.MakeVertex(    0,   0, 0 ) , Bubble1Points[0] ) 
vertical_line   = geompy.MakeLineTwoPnt( geompy.MakeVertex(    0,   0, 0 ) , Bubble1Points[-1] ) 


# Wire_1 = geompy.MakeWire([Curve_1, vertical_line, horizontal_line], 1e-07)
Bubble_to_be_cut = geompy.MakeFaceWires([ geompy.MakeWire([Curve_1, vertical_line, horizontal_line], 1e-07) ] , 1)

Domain_cut = geompy.MakeCutList( Domain_cut_uninflated, [Bubble_to_be_cut] )

geompy.addToStudy(Domain_cut,'Domain_cut')

# dR_refPlusEps = 1.2*dR_ref 
dR_refPlusEps = dR_ref 


Total_Scaled_Curve_Bubble1_Ref1_coords = []
Total_Scaled_Curve_Bubble1_Ref1_points = []
Total_Scaled_Curve_Bubble1_Ref1_coords = CurveCoordinatesScaling(Bubble1PointCoordinatesSorted, dR_refPlusEps)

Total_Scaled_Curve_Bubble1_Ref1_points = [geompy.MakeVertex(*p, 0) for p in Total_Scaled_Curve_Bubble1_Ref1_coords]
# Total_Scaled_Curve_Bubble1_Ref1_points = Total_Scaled_Curve_Bubble1_Ref1_points[:-1]


# i=0
# for item in Total_Scaled_Curve_Bubble1_Ref1_points:
#     if (Total_Scaled_Curve_Bubble1_Ref1_coords[i])[0]<0.04:
#         print(Total_Scaled_Curve_Bubble1_Ref1_coords[i])
#         Total_Scaled_Curve_Bubble1_Ref1_points.remove(item)
#     geompy.addToStudy(item,'point'+str(i))
#     i+=1

print(len(Total_Scaled_Curve_Bubble1_Ref1_coords))
print(len(Total_Scaled_Curve_Bubble1_Ref1_points))

for i in range(len(Total_Scaled_Curve_Bubble1_Ref1_points)):
    if (Total_Scaled_Curve_Bubble1_Ref1_coords[i])[0]<0.1:
        Total_Scaled_Curve_Bubble1_Ref1_points = Total_Scaled_Curve_Bubble1_Ref1_points[:i-1]


Total_Scaled_Curve_Bubble1_Ref1_points.append( geompy.MakeVertex( 0, Radius_needle, 0) )
Total_Scaled_Curve_Bubble1_Ref1 = geompy.MakeInterpol(Total_Scaled_Curve_Bubble1_Ref1_points, False)
geompy.addToStudy(Total_Scaled_Curve_Bubble1_Ref1,'Total_Scaled_Curve_Bubble1_Ref1')

pointOnSymmetry_coords = Bubble1PointCoordinatesSorted[0]

check_parameter = max(Bubble1PointCoordinatesSorted, key=lambda tup: tup[1])

# horizontal_line_square_around = geompy.MakeLineTwoPnt( geompy.MakeVertex( -Height_needle, Radius_needle + dR_refPlusEps + dH_ref, 0), geompy.MakeVertex( pointOnSymmetry_coords[0] + dR_refPlusEps + dH_ref, Radius_needle + dR_refPlusEps + dH_ref, 0) )
horizontal_line_square_around = geompy.MakeLineTwoPnt( geompy.MakeVertex( -Height_needle_modified, check_parameter[1] + dR_refPlusEps + dH_ref, 0), geompy.MakeVertex( pointOnSymmetry_coords[0] + dR_refPlusEps + dH_ref, check_parameter[1] + dR_refPlusEps + dH_ref, 0) )

# vertical_line_square_around = geompy.MakeLineTwoPnt( geompy.MakeVertex( pointOnSymmetry_coords[0] + dR_refPlusEps + dH_ref, Radius_needle + dR_refPlusEps + dH_ref, 0) , geompy.MakeVertex( pointOnSymmetry_coords[0] + dR_refPlusEps + dH_ref, 0, 0) )
vertical_line_square_around = geompy.MakeLineTwoPnt( geompy.MakeVertex( pointOnSymmetry_coords[0] + dR_refPlusEps + dH_ref, check_parameter[1] + dR_refPlusEps + dH_ref, 0) , geompy.MakeVertex( pointOnSymmetry_coords[0] + dR_refPlusEps + dH_ref, 0, 0) )

geompy.addToStudy(horizontal_line_square_around,'horizontal_line_square_around')
geompy.addToStudy(vertical_line_square_around,'vertical_line_square_around')
 

partition_line_list = [ Total_Scaled_Curve_Bubble1_Ref1, horizontal_line_square_around, vertical_line_square_around ]
partition_tool = geompy.MakeFuseList( partition_line_list, True, True)

Partition_1 = geompy.MakePartition([Domain_cut], [partition_tool], [], [], geompy.ShapeType["FACE"], 0, [], 0)

geompy.addToStudy(Partition_1,'Partition_1')

# x_tilt of type list [a,b,c]
def returnIDofShape( x_tilt, TypeofShape ):
    idShape     		 =  0.
    Shape_Point 		 = 	geompy.MakeVertex(*x_tilt)
    sLineShape  		 = 	geompy.GetShapesNearPoint(Partition_1, Shape_Point, geompy.ShapeType[ TypeofShape ])
    idShape 			 =  geompy.GetSubShapeID(Partition_1,	sLineShape  )    	
    return idShape

Height_tank = Height_tank - Height_needle

idWall = []
Wall = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idWall.append( returnIDofShape( [ 0,     Radius_needle-Needle_thickness+h_s, 0] , "EDGE" ) )
idWall.append( returnIDofShape( [              -h_s,     Radius_needle, 0] , "EDGE" ) )
idWall.append( returnIDofShape( [ Height_tank,     h_s, 0] , "EDGE" ) )
idWall.append( returnIDofShape( [-Height_needle_modified+h_s,       Radius_tank, 0] , "EDGE" ) )
Wall_union = geompy.UnionIDs( Wall, idWall )


idSymmetry = []
Symmetry = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idSymmetry.append( returnIDofShape( [  pointOnSymmetry_coords[0] + h_s,     0, 0] , "EDGE" ) )
idSymmetry.append( returnIDofShape( [ (Total_Scaled_Curve_Bubble1_Ref1_coords[0])[0] + h_s, (Total_Scaled_Curve_Bubble1_Ref1_coords[0])[1], 0 ] , "EDGE" ) )
idSymmetry.append( returnIDofShape( [ Height_tank-h_s,     0, 0] , "EDGE" ) )
Symmetry_union = geompy.UnionIDs( Symmetry, idSymmetry )


idInflatedBubble = []
pointOnInfaltedBubble_coords = Bubble1PointCoordinatesSorted[5]
InflatedBubble = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idInflatedBubble.append( returnIDofShape( [ pointOnInfaltedBubble_coords[0],  pointOnInfaltedBubble_coords[1], 0], "EDGE" ) )
InflatedBubble_union = geompy.UnionIDs( InflatedBubble, idInflatedBubble )


idAmbient = []
Ambient = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idAmbient.append( returnIDofShape( [    -Height_needle_modified, Radius_needle+h_s, 0] , "EDGE" ) )
idAmbient.append( returnIDofShape( [    -Height_needle_modified,   Radius_tank-h_s, 0] , "EDGE" ) )
Ambient_union = geompy.UnionIDs( Ambient, idAmbient )

# geompy.addToStudy(Ambient, "Ambient")
# geompy.addToStudy(InflatedBubble, "InflatedBubble")
# geompy.addToStudy(Wall, "Wall")
# geompy.addToStudy(Symmetry, "Symmetry")


idverticals1 = []
verticals1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idverticals1.append( idWall[0] )
# idverticals1.append( idSymmetry[0] )
verticals1_union = geompy.UnionIDs( verticals1, idverticals1 )

idhorizontals1 = []
pointOnScaledCurve_coords = Total_Scaled_Curve_Bubble1_Ref1_coords[5]
horizontals1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
idhorizontals1.append( idInflatedBubble[0] )
# idhorizontals1.append( returnIDofShape( [ pointOnScaledCurve_coords[0],  pointOnScaledCurve_coords[1], 0], "EDGE" ) )
horizontals1_union = geompy.UnionIDs( horizontals1, idhorizontals1 )

idquad_str = []
quad_str = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idquad_str.append( returnIDofShape( [  pointOnSymmetry_coords[0] + h_s,     h_s, 0] , "FACE" ) )
quad_str_union = geompy.UnionIDs( quad_str, idquad_str )

idquad_unstr = []
quad_unstr = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
idquad_unstr.append( returnIDofShape( [ (Total_Scaled_Curve_Bubble1_Ref1_coords[0])[0] + h_s, (Total_Scaled_Curve_Bubble1_Ref1_coords[0])[1]+h_s, 0 ] , "FACE" ) )
quad_str_union = geompy.UnionIDs( quad_unstr, idquad_unstr )



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
NETGEN_1D_2D_global   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
MeshParameters(NETGEN_1D_2D_global  ,Main_maxSize_element, Main_minSize_element, 0.3)

Mesh_1.Segment(geom=verticals1).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings,[])
Mesh_1.Segment(geom=verticals1).PropagationOfDistribution()

Mesh_1.Segment(geom=horizontals1).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings, [])
Mesh_1.Segment(geom=horizontals1).PropagationOfDistribution()

SM4 = Mesh_1.Quadrangle(geom = quad_str)
Mesh_1.Segment(geom=quad_str).NumberOfSegments(8)

SM6   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = quad_unstr)
MeshParameters(SM6  ,maxSquareElementLength, minSquareElementLength, 0.3)

# GSM0 = SM0.GetSubMesh()

# GSM2 = SM2.GetSubMesh()

GSM4 = SM4.GetSubMesh()

GSM6 = SM6.GetSubMesh()

# isDone = Mesh_1.SetMeshOrder( [   [ GSM0, GSM2, GSM4, GSM6] ])
isDone = Mesh_1.SetMeshOrder( [   [ GSM4, GSM6] ])

Wall_1                      =  Mesh_1.GroupOnGeom( Wall     ,'Wall'      ,SMESH.EDGE )
Symmetry_1                  =  Mesh_1.GroupOnGeom( Symmetry ,'Symmetry'  ,SMESH.EDGE )
InflatedBubble_1            =  Mesh_1.GroupOnGeom( InflatedBubble   ,'InflatedBubble'    ,SMESH.EDGE )
Ambient_1                   =  Mesh_1.GroupOnGeom( Ambient  ,'Ambient'   ,SMESH.EDGE )

isDone = Mesh_1.Compute()



# it likes cut 1 more since the cutting lines lie vertically to the moving surface when SECOND if condition is fulfilled
# it likes cut 0 more since the cutting lines lie vertically to the moving surface when FIRST if condition is fulfilled

isDone = Mesh_1.SplitQuadObject( Mesh_1, 0)
# FUTURE impovement: i can manually tell the program to cut the bubble wall quadrangle with cut 0 
# and the rest with cut 1



# isDone = Mesh_1.QuadTo4Tri( )
try:
  Mesh_1.ExportUNV( r'./UnBounded.unv' )
  pass
except:
  print('ExportUNV() failed. Invalid file name?')


#     Mesh_1.Segment(geom=verticals1).PropagationOfDistribution()



#     Mesh_1.Segment(geom=horizontals1).StartEndLength (startLengthOnSurroundings,endLengthOnSurroundings, idhorizontals1)
#     Mesh_1.Segment(geom=horizontals1).PropagationOfDistribution()

#     Mesh_1.Segment(geom=horizontals2).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings,[])
#     Mesh_1.Segment(geom=horizontals2).PropagationOfDistribution()


#     SM4 = Mesh_1.Quadrangle(geom = quad1)
#     Mesh_1.Segment(geom=quad1).NumberOfSegments(8)

#     SM5 = Mesh_1.Quadrangle(geom = quad2)
#     Mesh_1.Segment(geom=quad2).NumberOfSegments(8)


#     NETGEN_1D_2D_global   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
#     SM6   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = square)
#     MeshParameters(NETGEN_1D_2D_global  ,Main_maxSize_element, Main_minSize_element, 0.3)
#     MeshParameters(SM6  ,maxSquareElementLength, minSquareElementLength, 0.3)


#     GSM4 = SM4.GetSubMesh()
#     GSM5 = SM5.GetSubMesh()
#     GSM6 = SM6.GetSubMesh()

#     isDone = Mesh_1.SetMeshOrder( [   [ GSM4,
#                                         GSM5,
#                                         GSM6] ])

#     Wall_1                      =  Mesh_1.GroupOnGeom( Wall     ,'Wall'      ,SMESH.EDGE )
#     Symmetry_1                  =  Mesh_1.GroupOnGeom( Symmetry ,'Symmetry'  ,SMESH.EDGE )
#     InflatedBubble_1            =  Mesh_1.GroupOnGeom( InflatedBubble   ,'InflatedBubble'    ,SMESH.EDGE )
#     Ambient_1                   =  Mesh_1.GroupOnGeom( Ambient  ,'Ambient'   ,SMESH.EDGE )
#     InternalEquidistributionR_1 =  Mesh_1.GroupOnGeom( InternalEquidistributionR    ,'InternalEq_R'      ,SMESH.EDGE )


# else:
#     removal_list = []
#     for item in Total_Scaled_Curve_Bubble1_Ref1_coords:
#         if item[0]<0 and item[1]<=1 :
#             removal_list.append(item)



#     for item in removal_list:
#         Total_Scaled_Curve_Bubble1_Ref1_coords.remove(item)


#     delta_dist = abs( (Total_Scaled_Curve_Bubble1_Ref1_coords[-1])[0] - (Total_Scaled_Curve_Bubble1_Ref1_coords[-2])[0] )

#     last_point_coords = []
#     last_point_coords = [(Total_Scaled_Curve_Bubble1_Ref1_coords[-1])[0]+delta_dist, 1]
#     Total_Scaled_Curve_Bubble1_Ref1_coords.append(last_point_coords)


#     # Total_Scaled_Curve_Bubble1_Ref1_coords.append([ dR_refPlusEps, 0 ])

#     Total_Scaled_Curve_Bubble1_Ref1_points = [geompy.MakeVertex(*p, 0) for p in Total_Scaled_Curve_Bubble1_Ref1_coords]
#     Total_Scaled_Curve_Bubble1_Ref1 = geompy.MakeInterpol(Total_Scaled_Curve_Bubble1_Ref1_points, False)
    
#     i=0
#     for item in Total_Scaled_Curve_Bubble1_Ref1_points:
#         geompy.addToStudy(item,str(i))
#         i+=1

#     partition_line_list = [ Total_Scaled_Curve_Bubble1_Ref1, horizontal_line_square_around, vertical_line_square_around ]


#     partition_tool = geompy.MakeFuseList( partition_line_list, True, True)

#     geompy.addToStudy(partition_tool,'partition_tool')
#     Partition_1 = geompy.MakePartition([Domain_cut], [partition_tool], [], [], geompy.ShapeType["FACE"], 0, [], 0)

#     Height_tank = Height_tank - Height_needle

#     geompy.addToStudy(Partition_1,'Partition_1')



#     # idWall = []
#     # Wall = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#     # idWall.append( returnIDofShape( [              -h_s,     Radius_needle, 0] , "EDGE" ) )
#     # idWall.append( returnIDofShape( [    -Height_needle+h_s, Radius_needle, 0] , "EDGE" ) )
#     # idWall.append( returnIDofShape( [ Height_tank,     h_s, 0] , "EDGE" ) )
#     # # idWall.append( returnIDofShape( [    -Height_needle, Radius_needle+h_s, 0] , "EDGE" ) )
#     # # idWall.append( returnIDofShape( [-Height_needle    ,   Radius_tank-h_s, 0] , "EDGE" ) )
#     # idWall.append( returnIDofShape( [-Height_needle+h_s,       Radius_tank, 0] , "EDGE" ) )
#     # Wall_union = geompy.UnionIDs( Wall, idWall )


#     # idSymmetry = []
#     # Symmetry = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#     # idSymmetry.append( returnIDofShape( [  pointOnSymmetry_coords[0] + h_s,     0, 0] , "EDGE" ) )
#     # idSymmetry.append( returnIDofShape( [ (Total_Scaled_Curve_Bubble1_Ref1_coords[0])[0] + h_s, (Total_Scaled_Curve_Bubble1_Ref1_coords[0])[1], 0 ] , "EDGE" ) )
#     # idSymmetry.append( returnIDofShape( [ Height_tank-h_s,     0, 0] , "EDGE" ) )
#     # Symmetry_union = geompy.UnionIDs( Symmetry, idSymmetry )


#     # pointOnInfaltedBubble_coords = Bubble1PointCoordinatesSorted[5]
#     # idInflatedBubble = []
#     # InflatedBubble = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#     # idInflatedBubble.append( returnIDofShape( [ pointOnInfaltedBubble_coords[0],  pointOnInfaltedBubble_coords[1], 0], "EDGE" ) )
#     # InflatedBubble_union = geompy.UnionIDs( InflatedBubble, idInflatedBubble )

#     # idAmbient = []
#     # Ambient = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#     # idAmbient.append( returnIDofShape( [    -Height_needle_modified, Radius_needle+h_s, 0] , "EDGE" ) )
#     # idAmbient.append( returnIDofShape( [-Height_needle_modified    ,   Radius_tank-h_s, 0] , "EDGE" ) )
#     # # idAmbient.append( returnIDofShape( [ Height_tank,     h_s, 0] , "EDGE" ) )
#     # Ambient_union = geompy.UnionIDs( Ambient, idAmbient )


#     # idhorizontals1 = []
#     # horizontals1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#     # idhorizontals1.append( idWall[0] )
#     # idhorizontals1.append( idSymmetry[0] )
#     # horizontals1_union = geompy.UnionIDs( horizontals1, idhorizontals1 )

#     # idAcrossTheBubble = []
#     # AcrossTheBubble = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#     # idAcrossTheBubble.append( returnIDofShape( [ (Total_Scaled_Curve_Bubble1_Ref1_coords[5])[0] , (Total_Scaled_Curve_Bubble1_Ref1_coords[5])[1], 0 ] , "EDGE" ) )
#     # AcrossTheBubble_union = geompy.UnionIDs( AcrossTheBubble, idAcrossTheBubble )


#     # idquad1 = []
#     # quad1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
#     # idquad1.append( returnIDofShape( [ (Total_Scaled_Curve_Bubble1_Ref1_coords[5])[0]-h_s, (Total_Scaled_Curve_Bubble1_Ref1_coords[5])[1]-h_s, 0] , "FACE" ) )
#     # quad1_union = geompy.UnionIDs( quad1, idquad1 )


#     # idsquare = []
#     # square = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
#     # idsquare.append( returnIDofShape( [    -Height_needle + h_s, Radius_needle+dR_refPlusEps+h_s, 0] , "FACE" ) )
#     # square_union = geompy.UnionIDs( square, idsquare )




#     import  SMESH, SALOMEDS
#     from salome.smesh import smeshBuilder

#     smesh = smeshBuilder.New()


#     def MeshParameters( Properties, maxsize: float, minsize: float, growthRate: float): 

#         params = Properties.Parameters()
#         params.SetMaxSize            (  maxsize  )
#         params.SetMinSize            (  minsize  )
#         params.SetSecondOrder        (   0 )
#         params.SetOptimize           (   1 )
#         params.SetFineness           (   2 )
#         params.SetGrowthRate         (  growthRate )
#         params.SetChordalError       (  -1 )
#         params.SetChordalErrorEnabled(   0 )
#         params.SetUseSurfaceCurvature(   1 )
#         params.SetFuseEdges          (   1 )
#         params.SetWorstElemMeasure   (   0 )
#         params.SetUseDelauney        ( 108 )
#         params.SetQuadAllowed        (   0 )
#         params.SetCheckChartBoundary (   0 )




#     Mesh_1 = smesh.Mesh(Partition_1)


#     Mesh_1.Segment(geom=AcrossTheBubble).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings,[])
#     Mesh_1.Segment(geom=AcrossTheBubble).PropagationOfDistribution()


#     Mesh_1.Segment(geom=horizontals1).StartEndLength (startLengthOnSurroundings,startLengthOnSurroundings, [] )
#     # Mesh_1.Segment(geom=horizontals1).PropagationOfDistribution()


#     SM4 = Mesh_1.Quadrangle(geom = quad1)
#     Mesh_1.Segment(geom=quad1).NumberOfSegments(8)


#     NETGEN_1D_2D_global   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D)
#     SM6   = Mesh_1.Triangle(algo = smeshBuilder.NETGEN_1D2D, geom = square)
#     MeshParameters(NETGEN_1D_2D_global  ,Main_maxSize_element, Main_minSize_element, 0.3)
#     MeshParameters(SM6  ,maxSquareElementLength, minSquareElementLength, 0.3)
    
#     GSM4 = SM4.GetSubMesh()
#     GSM6 = SM6.GetSubMesh()

#     isDone = Mesh_1.SetMeshOrder( [   [ GSM4, GSM6] ])



#     Wall_1           			=  Mesh_1.GroupOnGeom( Wall 	,'Wall' 	 ,SMESH.EDGE )
#     Symmetry_1       			=  Mesh_1.GroupOnGeom( Symmetry ,'Symmetry'	 ,SMESH.EDGE )
#     InflatedBubble_1 			=  Mesh_1.GroupOnGeom( InflatedBubble 	,'InflatedBubble' 	 ,SMESH.EDGE )
#     Ambient_1        			=  Mesh_1.GroupOnGeom( Ambient 	,'Ambient' 	 ,SMESH.EDGE )


# isDone = Mesh_1.Compute()



# # it likes cut 1 more since the cutting lines lie vertically to the moving surface when SECOND if condition is fulfilled
# # it likes cut 0 more since the cutting lines lie vertically to the moving surface when FIRST if condition is fulfilled

# isDone = Mesh_1.SplitQuadObject( Mesh_1, 0)
# # FUTURE impovement: i can manually tell the program to cut the bubble wall quadrangle with cut 0 
# # and the rest with cut 1



# # isDone = Mesh_1.QuadTo4Tri( )
# try:
#   Mesh_1.ExportUNV( r'./UnBounded.unv' )
#   pass
# except:
#   print('ExportUNV() failed. Invalid file name?')



# if salome.sg.hasDesktop():
#   salome.sg.updateObjBrowser()


# import os

# try:
#     from killSalomeWithPort import killMyPort
#     killMyPort(os.getenv('NSPORT'))
# except:
#     pass