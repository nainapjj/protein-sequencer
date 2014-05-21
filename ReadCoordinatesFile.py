import ast

def readCoordinateFile(fileName):
    with open(fileName) as f:
        coordsByLine = f.readlines()
    
    coordinateSets = []
    
    for coordLine in coordsByLine:
        coordinateSets.append(ast.literal_eval(coordLine))
    
    return coordinateSets