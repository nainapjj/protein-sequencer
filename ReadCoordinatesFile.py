import ast
import numpy

import StuartFindAverageDistance

def readCoordinateFileOld(fileName):
    with open(fileName, 'r') as f:
        coordsByLine = f.readlines()

    coordinateSets = []
    
    for coordLine in coordsByLine:
        coordinateSets.append(ast.literal_eval(coordLine))
    
    return coordinateSets

def readCoordinateFile(fileName):
    with open(fileName, 'r') as f:
        coordsByLine = f.readlines()
    
    usedAmino = ast.literal_eval(coordsByLine[0])
    coordsByLine = coordsByLine[1]
    coordinateSets = ast.literal_eval(coordsByLine)
    
    return {"usedAmino": usedAmino, "coordSet": coordinateSets}

def findAverageLengthForEachSet(coordinateSets, usedAmino):
    distances = []
    
    for coordinateSet in coordinateSets:
        distances.append(StuartFindAverageDistance.LengthSet(coordinateSet, usedAmino))
    
    return distances

def findLengthForCoordinateSet(coordinateSet, usedAmino):
    # Attempt to find the Euclidean distance between each
    # subsequent coordinates
    distances = [-1] * (len(coordinateSet) - 1)

    for index in range(len(coordinateSet) - 1):
        if usedAmino[index+1] - usedAmino[index] == 1:
            coord_a = numpy.array(coordinateSet[index])
            coord_b = numpy.array(coordinateSet[index+1])

        distances[index] = numpy.linalg.norm(coord_a - coord_b)

    return distances
