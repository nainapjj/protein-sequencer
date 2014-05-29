import ast
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