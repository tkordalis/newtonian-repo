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



blockage_ratio = 0.01

# Domain Construction

RSphere1 		= 1.0
RSphere2 		= 1.0

Radius_tank  	= RSphere1/blockage_ratio
Height_tank 	= 200

Sphere_position = 0.


# NumSegmentsOnSphere = 200
NumSegmentsOnSphere = 40

Element_size_on_Sphere = 3.1415926535/NumSegmentsOnSphere



dR_ref1 = 3*Element_size_on_Sphere
dR_ref2 = 2*dR_ref1
R_refinement1_Sphere1 = RSphere1 + dR_ref1
R_refinement1_Sphere2 = RSphere2 + dR_ref1
R_refinement2_Sphere1 = RSphere1 + dR_ref2
R_refinement2_Sphere2 = RSphere2 + dR_ref2


Sphere_position + RSphere2 + dR_ref2 

# if (2*Sphere_position) <= R_refinement2_Sphere1 + R_refinement2_Sphere2:
# 	print("Sphere REFINEMENT ZONES are Touching")



ellipse_Minor_Radius = dR_ref2 + 1
ellipse_Major_Radius = dR_ref2 + 1

ellipse_position 	 = 0.



outer_ellipse_Major_Radius = ellipse_Major_Radius + 1.0
outer_ellipse_Minor_Radius = ellipse_Minor_Radius + 1.0


h_s = 0.0001


Main_maxSize_element = 4
Main_minSize_element = 0.1

Netgen_Params = []

check_drRef1_compatibility = int(dR_ref1/Element_size_on_Sphere)

if check_drRef1_compatibility<3:
	print('Segment length from sphere does not fit the symmetry egde')
	print(check_drRef1_compatibility)
	sys.exit()

# Fine Mesh parameters
Netgen_Params.append([ 2.*Element_size_on_Sphere,  2.*Element_size_on_Sphere, 0.5])

Netgen_Params.append([ 3*Element_size_on_Sphere,  3*Element_size_on_Sphere, 0.1])

Netgen_Params.append([ 4*Element_size_on_Sphere,  4*Element_size_on_Sphere, 0.1])
Netgen_Params.append([ 6*Element_size_on_Sphere,  6*Element_size_on_Sphere, 0.1])


# # --------------------------- End of Geometry --------------------------- #

