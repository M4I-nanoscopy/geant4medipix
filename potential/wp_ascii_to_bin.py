# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:35:15 2014

@author: armin
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
# number of grid points for different sizes
# 55um 300
# 110um 300 100x100x100 xy=4 z=3
# 110um 500 100x100x200   xy=4 z=2.5
# 110um 1000 100x100x400 xy=4 z=2.5
#
#
SENSOR=110
THICKN=1000
PATH = '/home/davkra/analysis/comsol/'

#FILES  = ["",
#	  "3d_regular_grid_weighting_field220_300_fine.txt",
#	  "",
#	  "3d_regular_grid_weighting_field110_300_fine.txt",
#          "3d_regular_grid_weighting_field110_500_fine.txt",
#	    "3d_regular_grid_weighting_field110_1000_fine.txt",
#          "3d_regular_grid_weighting_field55_200_fine.txt",
#          "3d_regular_grid_weighting_field55_300_fine.txt",
#          "3d_regular_grid_weighting_field55_500_fine.txt",
#          "3d_regular_grid_weighting_field55_1000_fine.txt"]

start = "3d_regular_grid_weighting_field"
end   = "_fine.txt"
# select parameters for pixel pitch
if SENSOR == 220:
    gridSizeFineXY=2.5
    nbGridPointsXY = 200
    files = start+str(SENSOR)
elif SENSOR == 110:
    gridSizeFineXY=2.5
    nbGridPointsXY = 200
    files = start+str(SENSOR)

elif SENSOR == 55:
    gridSizeFineXY = 2
    nbGridPointsXY = 200
    files = start+str(SENSOR)
else:
    sys.exit('Only 55, 110 or 220 are allowed pixel pitches!')

# select parameters for sensor thickness
if THICKN == 300:
    nbGridPointsZ = 150
    gridSizeFineZ = 2
    fname = files+'_'+str(THICKN)+end
elif THICKN == 500:
    nbGridPointsXY = 100
    nbGridPointsZ = 200
    gridSizeFineZ = 2.5
    fname = files+'_'+str(THICKN)+end
elif THICKN == 1000:
    nbGridPointsXY = 100
    gridSizeFineXY = 4
    nbGridPointsZ = 500
    gridSizeFineZ = 2.
    fname = files+'_'+str(THICKN)+end
elif THICKN == 200:
    nbGridPointsZ = 100
    gridSizeFineZ = 2.
    fname = files+'_'+str(THICKN)+end
else:
    sys.exit('Only 200, 300, 500 and 1000 are allowed thicknesses!')

halfgridz = gridSizeFineZ/2.

pot = np.zeros(shape=(nbGridPointsXY**2*nbGridPointsZ+4))

# read ascii file
print "start reading file..." + fname
f = open(os.path.join(PATH,fname), 'r')

#f = np.loadtxt('3d_regular_grid_weighting_field110_fine.txt', comments='%')


#store grid parameters in the file
pot[0] = gridSizeFineXY
pot[1] = gridSizeFineZ
pot[2] = nbGridPointsXY
pot[3] = nbGridPointsZ

print pot[0:4]

testarray = []
# get rid of the first 9 lines to skip checking for comments
for x in xrange(9):
    f.readline()

offset = nbGridPointsXY/2 - 0.5

for line in f:
    #if not line.startswith('%'):
    data = [splits for splits in line.split() if splits is not ""]
    #x = line.split('\t')
    x = int(offset+(np.float64(data[0])/gridSizeFineXY))
    y = int(offset+(np.float64(data[1])/gridSizeFineXY))
    z = int((np.float64(data[2])-halfgridz)/gridSizeFineZ)
    d = np.float64(data[3])
    #testarray.append((x+nbGridPointsXY*y+nbGridPointsXY*nbGridPointsXY*z)+4)
    pot[(x+nbGridPointsXY*y+nbGridPointsXY**2*z)+4] = d

# swap txt with bin in file name
base = os.path.splitext(fname)[0]
fname = base + ".bin"

# write file
print "Start writing file..." + fname
pot.tofile(os.path.join(fname))
print "...done"
