#!/usr/bin/env python

#> \file 
#> \author David Ladd, Reused: Hashem Yousefi 
#> \brief This is an example to use linear fitting to fit the beginning of linear heart tube.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): 
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>
### to be used for segmenting embryonic heart and fitting with an initial meshes.
#<

import sys, os
import exfile
import numpy
from numpy import linalg
import math
import random

# Intialise OpenCMISS/iron 
from opencmiss.iron import iron

# defining the output file to be written in the ExDataFile
def writeExdataFile(filename,dataPointLocations,dataErrorVector,dataErrorDistance,offset):
    "Writes data points to an exdata file"

    numberOfDimensions = dataPointLocations[1].shape[0]
    try:
        f = open(filename,"w")
        if numberOfDimensions == 1:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=2, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=3, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 2:
            header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=3, #Derivatives=0, #Versions=1
  y.  Value index=4, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=5, #Derivatives=0, #Versions=1
'''
        elif numberOfDimensions == 3:
             header = '''Group name: DataPoints
 #Fields=3
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=1, #Derivatives=0, #Versions=1
  y.  Value index=2, #Derivatives=0, #Versions=1
  x.  Value index=3, #Derivatives=0, #Versions=1
 2) data_error, field, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  x.  Value index=4, #Derivatives=0, #Versions=1
  y.  Value index=5, #Derivatives=0, #Versions=1
  z.  Value index=6, #Derivatives=0, #Versions=1
 3) data_distance, field, real, #Components=1
  1.  Value index=7, #Derivatives=0, #Versions=1
'''
        f.write(header)

        numberOfDataPoints = len(dataPointLocations)
        for i in range(numberOfDataPoints):
            line = " Node: " + str(offset+i+1) + '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataPointLocations[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataErrorVector[i,j]) + '\t'
                f.write(line)
            line = '\n'
            f.write(line)
            line = ' ' + str(dataErrorDistance[i])
            f.write(line)
            line = '\n'
            f.write(line)
        f.close()
            
    except IOError:
        print ('Could not open file: ' + filename)

#=================================================================
# Control Panel
#=================================================================
# set the number of elements and the number of nodes for the cylinder 
numberOfDimensions = 3
numberOfGaussXi = 3 
numberOfCircumfrentialElementsPerQuarter = 2
numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfLengthElements = 8
numberOfLengthNodes = numberOfLengthElements+1
numberOfWallElements = 1
numberOfWallNodes = numberOfWallElements+1
origin = [0.0,0.0,0.0]
meshOrigin = [0.0,0.0,0.0]
print "mesh resolution and parameters fixed"

# The number of data points which are digitised from the heart segments 
# fix interior nodes so that fitting only applies to surface
# If start iteration > 1, read in geometry from a previous fit iteration
numberOfDataPoints = 3
numberOfIterations = 1
fixInterior = True
zeroTolerance = 0.00001
hermite = True
# Set Sobolev smoothing parameters
tau = 0.00005
kappa = 0.00001
iteration = 1
if iteration > 1:
    exfileMesh = True
    exnode = exfile.Exnode("DeformedGeometry" + str(iteration-1) + ".part0.exnode")
    exelem = exfile.Exelem("UndeformedGeometry.part0.exelem")
else:
    exfileMesh = False
print "other parameters setted up "

#=================================================================
# 
#=================================================================
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    dataPointFieldUserNumber,
    materialFieldUserNumber,
    analyticFieldUserNumber,
    dependentDataFieldUserNumber,
    dataPointsUserNumber,
    dataProjectionUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,19)

# Get the computational nodes information
print dir(iron),'\n\n'
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC CS
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# define a basis 
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfDimensions
if hermite:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
else:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()
print "CS, Region and basis setted up"

#=================================================================
# Mesh
#=================================================================
# creating the number of elements and the mesh origins ... and/or
# Start the creation of a manually generated mesh in the region
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements

print "numberOfElements = ", numberOfElements
print "numberOfNodes = ", numberOfNodes


if (exfileMesh):
    # Read previous mesh
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, numberOfDimensions)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(exelem.num_elements)
    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, exnode.num_nodes)
    nodes.CreateFinish()
    # Define elements for the mesh
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    for elem in exelem.elements:
        elements.NodesSet(elem.number, elem.nodes)
    elements.CreateFinish()
    mesh.CreateFinish()
else:
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber,region,3)
    mesh.origin = meshOrigin
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(numberOfElements)
# Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region,numberOfNodes)
    nodes.CreateFinish()
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    elementNumber = 0
    for wallElementIdx in range(1,numberOfWallElements+1):
       for lengthElementIdx in range(1,numberOfLengthElements+1):
            for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
                elementNumber = elementNumber + 1
                localNode1 = circumfrentialElementIdx + (lengthElementIdx - 1)*numberOfCircumfrentialElements + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                if circumfrentialElementIdx == numberOfCircumfrentialElements:
                    localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                        (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
                else: 
                    localNode2 = localNode1 + 1
                localNode3 = localNode1 + numberOfCircumfrentialNodes
                localNode4 = localNode2 + numberOfCircumfrentialNodes
                localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
                localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
#		print "Element Number = ",elementNumber
#         	print "Node numbers of the element", localNode1, localNode2, localNode3, localNode4, localNode5, localNode6, localNode7, localNode8 
                elements.NodesSet(elementNumber,localNodes)  
    elements.CreateFinish()
    mesh.CreateFinish() 


# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

print "mesh decomposition finished"

#=================================================================
# Geometric Field
#=================================================================
# the location of  nodes for the mesh  
manualNodePoints = numpy.zeros((9,8,3,2))

manualNodePoints[0,0,:,0] = [1190,685,-25]
manualNodePoints[0,1,:,0] = [1090,730,-25]
manualNodePoints[0,2,:,0] = [920,745,-25]
manualNodePoints[0,3,:,0] = [795,715,-25]
manualNodePoints[0,4,:,0] = [745,625,-25]
manualNodePoints[0,5,:,0] = [700,550,-25]
manualNodePoints[0,6,:,0] = [840,535,-25]
manualNodePoints[0,7,:,0] = [1010,585,-25]

manualNodePoints[1,0,:,0] = [1615,725,375]
manualNodePoints[1,1,:,0] = [1560,775,375]
manualNodePoints[1,2,:,0] = [1430,820,375]
manualNodePoints[1,3,:,0] = [1310,788,375]
manualNodePoints[1,4,:,0] = [1270,735,375]
manualNodePoints[1,5,:,0] = [1375,675,375]
manualNodePoints[1,6,:,0] = [1475,650,375]
manualNodePoints[1,7,:,0] = [1590,680,375]

manualNodePoints[2,7,:,0] = [1466,558,710]
manualNodePoints[2,6,:,0] = [1401,539,681]
manualNodePoints[2,5,:,0] = [1342,538,656]
manualNodePoints[2,4,:,0] = [1281,563,630]
manualNodePoints[2,3,:,0] = [1309,625,661]
manualNodePoints[2,2,:,0] = [1373,673,684]
manualNodePoints[2,1,:,0] = [1438,708,702]
manualNodePoints[2,0,:,0] = [1524,671,736]

'''
manualNodePoints[2,7,:,0] = [1449,575,760]
manualNodePoints[2,6,:,0] = [1370,550,720]
manualNodePoints[2,5,:,0] = [1325,545,680]
manualNodePoints[2,4,:,0] = [1253,561,637]
manualNodePoints[2,3,:,0] = [1295,630,680]
manualNodePoints[2,2,:,0] = [1360,685,720]
manualNodePoints[2,1,:,0] = [1430,730,760]
manualNodePoints[2,0,:,0] = [1520,700,810]
'''

manualNodePoints[3,7,:,0] = [1419,461,1039]
manualNodePoints[3,6,:,0] = [1369,375,902]
manualNodePoints[3,5,:,0] = [1322,409,804]
manualNodePoints[3,4,:,0] = [1261,470,722]
manualNodePoints[3,3,:,0] = [1207,563,744]
manualNodePoints[3,2,:,0] = [1197,579,853]
manualNodePoints[3,1,:,0] = [1245,558,931]
manualNodePoints[3,0,:,0] = [1312,557,1040]

manualNodePoints[4,7,:,0] = [1113,306,1051]
manualNodePoints[4,6,:,0] = [1217,354,925]
manualNodePoints[4,5,:,0] = [1181,380,825]
manualNodePoints[4,4,:,0] = [1143,415,760]
manualNodePoints[4,3,:,0] = [1074,557,749]
manualNodePoints[4,2,:,0] = [1024,573,829]
manualNodePoints[4,1,:,0] = [1043,514,950]
manualNodePoints[4,0,:,0] = [1025,421,1013]

manualNodePoints[5,6,:,0] = [780,255,672]
manualNodePoints[5,5,:,0] = [902,245,500]
manualNodePoints[5,4,:,0] = [917,290,400]
manualNodePoints[5,3,:,0] = [952,438,560]
manualNodePoints[5,2,:,0] = [865,500,700]
manualNodePoints[5,1,:,0] = [850,570,784]
manualNodePoints[5,0,:,0] = [710,435,900]
manualNodePoints[5,7,:,0] = [766,223,840]

manualNodePoints[6,0,:,0] = [496,526,555]
manualNodePoints[6,1,:,0] = [595,605,600]
manualNodePoints[6,2,:,0] = [700,632,575]
manualNodePoints[6,3,:,0] = [588,592,512]
manualNodePoints[6,4,:,0] = [510,544,470]
manualNodePoints[6,5,:,0] = [357,492,435]
manualNodePoints[6,6,:,0] = [256,412,431]
manualNodePoints[6,7,:,0] = [385,463,517]

manualNodePoints[7,1,:,0] = [570,850,855]
manualNodePoints[7,2,:,0] = [515,900,855]
manualNodePoints[7,3,:,0] = [350,885,855]
manualNodePoints[7,4,:,0] = [160,800,855]
manualNodePoints[7,5,:,0] = [48.0,700,855]
manualNodePoints[7,6,:,0] = [110,670,855]
manualNodePoints[7,7,:,0] = [270,700,855]
manualNodePoints[7,0,:,0] = [450,775,855]

manualNodePoints[8,1,:,0] = [1030,560,1200]
manualNodePoints[8,2,:,0] = [770,610,1200]
manualNodePoints[8,3,:,0] = [600,670,1200]
manualNodePoints[8,4,:,0] = [230,750,1200]
manualNodePoints[8,5,:,0] = [170,600,1200]
manualNodePoints[8,6,:,0] = [290,450,1200]
manualNodePoints[8,7,:,0] = [658,438,1200]
manualNodePoints[8,0,:,0] = [900,475,1200]

# node positions of the outer surface ... 

manualNodePoints[0,7,:,1] = [1180,310,-25]
manualNodePoints[0,6,:,1] = [865,365,-25]
manualNodePoints[0,5,:,1] = [660,535,-25]
manualNodePoints[0,4,:,1] = [650,650,-25]
manualNodePoints[0,3,:,1] = [740,760,-25]
manualNodePoints[0,2,:,1] = [860,810,-25]
manualNodePoints[0,1,:,1] = [1075,810,-25]
manualNodePoints[0,0,:,1] = [1300,700,-25]

manualNodePoints[1,7,:,1] = [1750,540,375]
manualNodePoints[1,6,:,1] = [1409,485,375]
manualNodePoints[1,5,:,1] = [1100,560,375]
manualNodePoints[1,4,:,1] = [813,701,368]
manualNodePoints[1,3,:,1] = [1000,900,375]
manualNodePoints[1,2,:,1] = [1365,880,375]
manualNodePoints[1,1,:,1] = [1670,863,375]
manualNodePoints[1,0,:,1] = [1880,710,375]

manualNodePoints[2,7,:,1] = [1720,485,710]
manualNodePoints[2,6,:,1] = [1572,407,640]
manualNodePoints[2,5,:,1] = [1357,540,571]
manualNodePoints[2,4,:,1] = [1160,568,561]
manualNodePoints[2,3,:,1] = [1079,697,604]
manualNodePoints[2,2,:,1] = [1332,787,718]
manualNodePoints[2,1,:,1] = [1541,846,710]
manualNodePoints[2,0,:,1] = [1825,616,862]

'''
manualNodePoints[2,7,:,1] = [1667,450,866]
manualNodePoints[2,6,:,1] = [1500,380,750]
manualNodePoints[2,5,:,1] = [1375,323,660]
manualNodePoints[2,4,:,1] = [1238,588,560]
manualNodePoints[2,3,:,1] = [1220,742,660]
manualNodePoints[2,2,:,1] = [1341,789,750]
manualNodePoints[2,1,:,1] = [1502,837,868]
manualNodePoints[2,0,:,1] = [1807,619,1074]
'''
manualNodePoints[3,7,:,1] = [1582,282,1071]
manualNodePoints[3,6,:,1] = [1413,211,772]
manualNodePoints[3,5,:,1] = [1294,181,475]
manualNodePoints[3,4,:,1] = [1091,522,546]
manualNodePoints[3,3,:,1] = [979,696,587]
manualNodePoints[3,2,:,1] = [1133,695,848]
manualNodePoints[3,1,:,1] = [1298,723,1070]
manualNodePoints[3,0,:,1] = [1728,642,1239]

manualNodePoints[4,7,:,1] = [1135,107,1058]
manualNodePoints[4,6,:,1] = [1020,105,704]
manualNodePoints[4,5,:,1] = [1055,166,464]
manualNodePoints[4,4,:,1] = [999,412,376]
manualNodePoints[4,3,:,1] = [892,678,533]
manualNodePoints[4,2,:,1] = [970,759,731]
manualNodePoints[4,1,:,1] = [1083,571,1048]
manualNodePoints[4,0,:,1] = [1261,270,1383]

manualNodePoints[5,7,:,1] = [739,102,1093]
manualNodePoints[5,6,:,1] = [782,83,674]
manualNodePoints[5,5,:,1] = [687,131,268]
manualNodePoints[5,4,:,1] = [764,275,113]
manualNodePoints[5,3,:,1] = [780,655,400]
manualNodePoints[5,2,:,1] = [853,785,606]
manualNodePoints[5,1,:,1] = [820,718,872]
manualNodePoints[5,0,:,1] = [574,445,993]

manualNodePoints[6,7,:,1] = [181,263,450]
manualNodePoints[6,6,:,1] = [-39,384,209]
manualNodePoints[6,5,:,1] = [264,693,257]
manualNodePoints[6,4,:,1] = [520,746,414]
manualNodePoints[6,3,:,1] = [653,804,507]
manualNodePoints[6,2,:,1] = [707,774,661]
manualNodePoints[6,1,:,1] = [608,586,714]
manualNodePoints[6,0,:,1] = [369,426,614]

manualNodePoints[7,0,:,1] = [545,667,855]
manualNodePoints[7,1,:,1] = [655,872,855]
manualNodePoints[7,2,:,1] = [525,940,855]
manualNodePoints[7,3,:,1] = [320,950,855]
manualNodePoints[7,4,:,1] = [75,840,855]
manualNodePoints[7,5,:,1] = [-70,710,855]
manualNodePoints[7,6,:,1] = [90,550,855]
manualNodePoints[7,7,:,1] = [242,590,855]

manualNodePoints[8,7,:,1] = [675,400,1200]
manualNodePoints[8,6,:,1] = [305,370,1200]
manualNodePoints[8,5,:,1] = [70,580,1200]
manualNodePoints[8,4,:,1] = [210,795,1200]
manualNodePoints[8,3,:,1] = [635,700,1200]
manualNodePoints[8,2,:,1] = [834,660,1200]
manualNodePoints[8,1,:,1] = [1150,620,1200]
manualNodePoints[8,0,:,1] = [940,430,1200]


#calculating the derivatives 
difference = numpy.zeros((9,8,3,2))
differenceAverage = numpy.zeros((9,8,3,2))
circumDeriv = numpy.zeros((9,8,3,2))
directDeriv = numpy.zeros((9,8,3,2))
lengthDeriv = numpy.zeros((9,8,3,2))
#circumferential derivative to be calculated 
for k in range (2):
    for j in range (9):
        for i in range (8):
            if (i<7):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i+1,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,0,m,k]-manualNodePoints[j,7,m,k]
for k in range (2):
    for j in range (9):
        for i in range (8):
            if (i<7):
                for m in range (3):
                    differenceAverage[j,i+1,m,k]=(difference[j,i+1,m,k]+difference[j,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,0,m,k]=(difference[j,0,m,k]+difference[j,7,m,k])/2
for k in range (2):
    for j in range (9):
        for i in range (8):
            for m in range (3):
                circumDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# derivative of the length direction
for k in range (2):
    for i in range (8):
        for j in range (9):
            if (j<8):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j+1,i,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i,m,k]-manualNodePoints[j-1,i,m,k]
for k in range (2):
    for i in range (8):
        for j in range (9):
            if (j == 0):
                for m in range (3): 
                    differenceAverage[j,i,m,k]=difference[j,i,m,k]
            if (j<8):
                for m in range (3):
                    differenceAverage[j+1,i,m,k]=(difference[j,i,m,k]+difference[j+1,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,i,m,k]=difference[j-1,i,m,k]
for k in range (2):
    for j in range (9):
        for i in range (8):
            for m in range (3):
                lengthDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# the derivatives of the wall direction is defined in the below lines ... 
for i in range (8):
    for j in range (9):
        for m in range (3):
            for k in range (2):
                difference[j,i,m,k] = manualNodePoints[j,i,m,1] - manualNodePoints[j,i,m,0]
for i in range (8):
    for j in range (9):
        for k in range (2):
            for m in range (3):
                directDeriv[j,i,m,k] = difference[j,i,m,k]/math.sqrt(math.pow(difference[j,i,0,k],2) + math.pow(difference[j,i,1,k],2) + math.pow(difference[j,i,2,k],2))


# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
for dimension in range(3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimension+1,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Create the geometric field
for wallNodeIdx in range(1,numberOfWallNodes+1):
    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
            x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            xtangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            ytangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            ztangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            xnormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            ynormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            znormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            zxnormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            zynormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            zznormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,1,nodeNumber,1,x)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,1,nodeNumber,2,y)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,1,nodeNumber,3,z)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,ztangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,zxnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,zynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,zznormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,znormal)


# Update the geometric field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes

# Get or calculate geometric parameters
if (exfileMesh):
    # Read the geometric field from the exnode file
    geometricField.ParameterSetUpdateStart(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    for node_num in range(1, exnode.num_nodes + 1):
        version = 1
        derivative = 1
        for component in range(1, numberOfDimensions + 1):
            component_name = ["x", "y", "z"][component - 1]
            value = exnode.node_value("Coordinate", component_name, node_num, derivative)
            geometricField.ParameterSetUpdateNode(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    version, derivative, node_num, component, value)
    geometricField.ParameterSetUpdateFinish(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
else:
    # Create undeformed geometry from the generated mesh
#    mesh.GeometricParametersCalculate(geometricField)
    # Export undeformed mesh geometry
    print("Writing undeformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("UndeformedGeometry","FORTRAN")
    fields.ElementsExport("UndeformedGeometry","FORTRAN")
    fields.Finalise()

#=================================================================
# Data Points
#=================================================================
# Create the data points
dataPoints = iron.DataPoints()
dataPoints.CreateStart(dataPointsUserNumber,region,numberOfDataPoints)
dataPointLocations = numpy.zeros((numberOfDataPoints,3))
print("Number of data points: " + str(numberOfDataPoints))
# reading from a text file containing the point clouds   
with open("EpicardialFullPoints.txt", "r") as ins:
	arrayOfInputData = []
	for line in ins:
		arrayOfInputData.append(line)
x = 0
y = 0
z = 0
for i in range (numberOfDataPoints):
	for j in range (5):
		sample = arrayOfInputData[i*5 + j]
		if (math.fmod(j,5) == 1):
			x = float (sample[12:25])				
		elif (math.fmod(j,5) == 2):
			y = float (sample[12:25])
		elif (math.fmod(j,5) == 3):
			z = float (sample[12:17])
		dataPointLocations[i,:] = [x,y,z]
# Set up data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.PositionSet(dataPointId,dataList)
dataPoints.CreateFinish()
 
'''
# write data points to exdata file for CMGUI
offset = 0
writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,offset)
'''

#=================================================================
# Data Projection on Geometric Field
#=================================================================
print("Projecting data points onto geometric field")
candidateElements = range(1,numberOfElements+1)
candidateFaceNormals = iron.ElementNormalXiDirections.PLUS_XI3*numpy.ones(numberOfElements,dtype=numpy.int32)
# Set up data projection
dataProjection = iron.DataProjection()
dataProjection.CreateStart(dataProjectionUserNumber,dataPoints,geometricField,iron.FieldVariableTypes.U)
#dataProjection.projectionType = iron.DataProjectionProjectionTypes.ALL_ELEMENTS
dataProjection.projectionType = iron.DataProjectionProjectionTypes.BOUNDARY_FACES
dataProjection.ProjectionCandidateFacesSet(candidateElements,candidateFaceNormals)
dataProjection.ProjectionDataCandidateFacesSet([1,2,3],[1,2],[iron.ElementNormalXiDirections.PLUS_XI3,iron.ElementNormalXiDirections.PLUS_XI3])
dataProjection.CreateFinish()

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition topology for data projection
decomposition.TopologyDataProjectionCalculate()

# Output data projection results
dataProjection.ResultAnalysisOutput("")

# Cancel some projections

dataProjection.ProjectionCancelByExitTags([iron.DataProjectionExitTags.EXIT_TAG_MAX_ITERATION,iron.DataProjectionExitTags.EXIT_TAG_NO_ELEMENT])
dataProjection.ResultAnalysisOutput("")

dataProjection.ProjectionCancelByDistance(iron.DataProjectionDistanceRelations.GREATER_EQUAL,250.0)
dataProjection.ResultAnalysisOutput("")

dataProjection.ProjectionCancelByDataPoints([1021,1024])
dataProjection.ResultAnalysisOutput("")

# Output the .exdata file.                                           
dataErrorVector = numpy.zeros((numberOfDataPoints,3))
dataErrorDistance = numpy.zeros(numberOfDataPoints)
for elementIdx in range(1,numberOfElements+1):
    numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementIdx)
    for dataPointIdx in range(1,numberOfProjectedDataPoints+1):
        dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementIdx,dataPointIdx)
        errorVector = dataProjection.ResultProjectionVectorGet(dataPointNumber,3)
        print("data point number = ",dataPointNumber,", errorVector = ",errorVector)
        dataErrorVector[dataPointNumber-1,0]=errorVector[0]
        dataErrorVector[dataPointNumber-1,1]=errorVector[1]
        dataErrorVector[dataPointNumber-1,2]=errorVector[2]
        errorDistance = dataProjection.ResultDistanceGet(dataPointNumber)
        dataErrorDistance[dataPointNumber-1]=errorDistance
 
# write data points to exdata file for CMGUI
offset = 0
writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,dataErrorVector,dataErrorDistance,offset)

print("Projection complete")
 
 
#=================================================================
# Equations Set
#=================================================================
# Create vector fitting equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.DATA_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.DATA_POINT_FITTING, 
 							 iron.EquationsSetFittingSmoothingTypes.SOBOLEV_VALUE]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

#=================================================================
# Dependent Field
#=================================================================
# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,numberOfDimensions)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions)
equationsSet.DependentCreateFinish()
# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

# Initialise dependent field to undeformed geometric field
for component in range (1,numberOfDimensions+1):
    geometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, component)

#=================================================================
# Independent Field
#=================================================================
# Create data point field (independent field, with vector values stored at the data points)
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U,"data point vector")
independentField.VariableLabelSet(iron.FieldVariableTypes.V,"data point weight")
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()
# Initialise data point vector field to 0
#independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
# Initialise data point weight field to 1
#independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)
# loop over each element's data points and set independent field values to data point locations on surface of the sphere
for element in range(numberOfElements):
    elementId = element + 1
    elementDomain = decomposition.ElementDomainGet(elementId)
    if (elementDomain == computationalNodeNumber):
        numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementId)
        for dataPoint in range(numberOfProjectedDataPoints):
            dataPointId = dataPoint + 1
            dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementId,dataPointId)
            dataList = dataPoints.ValuesGet(dataPointNumber,3)
            # set data point field values
            for component in range(numberOfDimensions):
                componentId = component + 1
                dataPointNumberIndex = dataPointNumber - 1
                value = dataList[component]
                independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementId,dataPointId,componentId,value)

#=================================================================
# Material Field
#=================================================================
# Create material field (Sobolev parameters)
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Smoothing Parameters")
equationsSet.MaterialsCreateFinish()
# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================
# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.FULL
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================
# Create fitting problem
problem = iron.Problem()
#problemSpecification = [iron.ProblemClasses.FITTING,
#                        iron.ProblemTypes.DATA_FITTING,
#                        iron.ProblemSubtypes.DATA_POINT_VECTOR_STATIC_FITTING]
problemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.STATIC_FITTING]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE # NONE / MATRIX
#solver.outputType = iron.SolverOutputTypes.MATRIX # NONE / MATRIX
solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.linearType = iron.LinearSolverTypes.DIRECT
#solver.LibraryTypeSet(iron.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
#solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
solver.linearIterativeAbsoluteTolerance = 1.0E-10
solver.linearIterativeRelativeTolerance = 1.0E-05
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
#solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

version = 1
meshComponent = decomposition.MeshComponentGet()
# Fix the interior nodes- use to only apply fit to surface nodes
if (fixInterior):
    # first find which nodes are non-surface nodes
    for nodeIdx in range(1,9):
        nodeDomain = decomposition.NodeDomainGet(nodeIdx,meshComponent)
        if (nodeDomain == computationalNodeNumber):
            for componentIdx in range(1,numberOfDimensions+1):
                numberOfDerivatives=meshNodes.NumberOfDerivativesGet(nodeIdx)
                for derivativeIdx in range(1,numberOfDerivatives+1):
                    value=geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,
                                                               version,derivativeIdx,nodeIdx,componentIdx)
                    #print "Setting boundary condition for component, node, derivative, value: ",componentIdx,nodeIdx,derivativeIdx,value
                    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,
                                               version,derivativeIdx,nodeIdx,componentIdx,
                                               iron.BoundaryConditionsTypes.FIXED,value)

solverEquations.BoundaryConditionsCreateFinish()


#=================================================================
# S o l v e    a n d    E x p o r t    D a t a
#=================================================================
derivativeVector=[0.0,0.0,0.0,0.0]
for iteration in range (1,numberOfIterations+1):
    # Solve the problem
    print("Solving fitting problem, iteration: " + str(iteration))
    problem.Solve()
    # Normalise derivatives
    for nodeIdx in range(1,numberOfNodes+1):
      for derivativeIdx in [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]:
          length=0.0
          for componentIdx in range(1,4):
              derivativeVector[componentIdx]=dependentField.ParameterSetGetNode(iron.FieldVariableTypes.U,
                                                                                iron.FieldParameterSetTypes.VALUES,
                                                                                1,derivativeIdx,nodeIdx,componentIdx)
              length=length + derivativeVector[componentIdx]*derivativeVector[componentIdx]
          length=math.sqrt(length)
          for componentIdx in range(1,4):
              value=derivativeVector[componentIdx]/length
              dependentField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    1,derivativeIdx,nodeIdx,componentIdx,value)
    # Copy dependent field to geometric 
    for componentIdx in range(1,numberOfDimensions+1):
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx,geometricField,
                                                                iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                componentIdx)
    # Reproject
    dataProjection.DataPointsProjectionEvaluate(iron.FieldParameterSetTypes.VALUES)
    rmsError=dataProjection.ResultRMSErrorGet()
    print("RMS error = "+ str(rmsError))
    # Export fields
    print("Writing deformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()
#-----------------------------------------------------------------

iron.Finalise()
