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
cavity_edge = 1

h_s = 5e-05

NumSegmentsOnCavityEdge = 120

