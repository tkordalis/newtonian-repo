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



blockage_ratio = 0.02

# Domain Construction

RSphere1 		= 1.0
RSphere2 		= 1.0

Radius_tank  	= RSphere1/blockage_ratio
Height_tank 	= 200

Sphere_position = 0.


# NumSegmentsOnSphere = 200
NumSegmentsOnSphere  = 600
NumSegmentsOnAmbient = 30

Element_size_on_Sphere = 3.1415926535/NumSegmentsOnSphere

Element_size_on_Ambient= Radius_tank/NumSegmentsOnAmbient
Element_size_on_Ambient_cb = 1*Element_size_on_Ambient


dR_ref1 = 4*Element_size_on_Sphere
dR_ref2 = 3*dR_ref1
R_refinement1_Sphere1 = RSphere1 + dR_ref1
R_refinement2_Sphere1 = RSphere1 + dR_ref2

dZ_refAmb = 3.0*Element_size_on_Ambient_cb
dZ_refAmb = 0.5*Height_tank - dZ_refAmb


# if (2*Sphere_position) <= R_refinement2_Sphere1 + R_refinement2_Sphere2:
# 	print("Sphere REFINEMENT ZONES are Touching")



ellipse_Minor_Radius = R_refinement2_Sphere1 + 1
ellipse_Major_Radius = R_refinement2_Sphere1 + 1.1

ellipse_position 	 = 0.



outer_ellipse_Minor_Radius = ellipse_Minor_Radius + 1.0
outer_ellipse_Major_Radius = ellipse_Major_Radius + 1.0


h_s = 0.0001


Main_maxSize_element = 2
Main_minSize_element = 0.1

NodeDensityFunction_Sym   = '(2*t-0.99)^6+0.05' 


Netgen_Params = []

check_drRef1_compatibility = int(dR_ref1/Element_size_on_Sphere)


# Fine Mesh parameters
Netgen_Params.append([ 2*Element_size_on_Sphere,  2*Element_size_on_Sphere, 0.1])

Netgen_Params.append([ 4*Element_size_on_Sphere,  4*Element_size_on_Sphere, 0.1])
Netgen_Params.append([ 8*Element_size_on_Sphere,  8*Element_size_on_Sphere, 0.1])

# # --------------------------- End of Geometry --------------------------- #

