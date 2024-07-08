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

# from the photos shown in Fig. 4 they mention that the inner diameter is 0.58 mm
# and from the scale bar we can deduce that the outer diameter is 0.8 mm
# which means the wall thickness is 2 x 0.11 mm -> the ratio is 0.11/0.8 = 0.1375 constant throughout the simulations



# Domain Construction

Height_needle = 5
# Radius_needle = 1.31
Radius_needle = 1.448

Needle_thickness = Radius_needle-1.0
Needle_thickness08 = 0.95*Needle_thickness

Height_tank = 20
Height_tank = Height_tank + Height_needle
Radius_tank = 10

dR_ref = Needle_thickness08
dH_ref = 4


h_s = 5e-05

NumSegmentsOnInflatedBubble = 5
Main_maxSize_element = 10	
Main_minSize_element = 0.8

edgeSizeOnInflatedBubble = (Radius_needle-Needle_thickness)/NumSegmentsOnInflatedBubble

NumSegmentsOnSymmetry = int(dR_ref/edgeSizeOnInflatedBubble)

startLengthOnSurroundings =   edgeSizeOnInflatedBubble
# endLengthOnSurroundings   =   edgeSizeOnInflatedBubble
endLengthOnSurroundings   = edgeSizeOnInflatedBubble

minSquareElementLength = 4*startLengthOnSurroundings
maxSquareElementLength = 4*startLengthOnSurroundings
# minSquareElementLength = 8*startLengthOnSurroundings
# maxSquareElementLength = 8*startLengthOnSurroundings
