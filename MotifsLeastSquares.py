import scipy
import scipy.linalg
import scipy.optimize
import random
import re

import Constants

findResidualsCount = 0

# Groups: 1 - X, 3 - Y, 5 - Z
DIFFERENCE_COORDINATES_PTN = r"(?<=\()([\-\.0-9]+)(, )([\-\.0-9]+)(, )([\-\.0-9]+)"

def generateAminoAcidsDifferenceList():
    listOfAminoDifferences = []
    with open(Constants.differenceInAminoFile) as f:
        for line in f.readlines():
            coord_match = re.search(DIFFERENCE_COORDINATES_PTN, line)
            x = float(coord_match.group(1))
            y = float(coord_match.group(3))
            z = float(coord_match.group(5))
            listOfAminoDifferences.append((x, y, z))
    return listOfAminoDifferences

def getCoordinatesByAminoAcidIndex(aminoIndices, coordinates, usedAmino):
    v0 = usedAmino.index(aminoIndices[0])
    v1 = usedAmino.index(aminoIndices[1])
    v2 = usedAmino.index(aminoIndices[2])
    v3 = usedAmino.index(aminoIndices[3])
    aa1Coordinates = (coordinates[v0 * 3 + 0], coordinates[v0 * 3 + 1],
                      coordinates[v0 * 3 + 2])
    aa2Coordinates = (coordinates[v1 * 3 + 0], coordinates[v1 * 3 + 1],
                      coordinates[v1 * 3 + 2])
    aa3Coordinates = (coordinates[v2 * 3 + 0], coordinates[v2 * 3 + 1],
                      coordinates[v2 * 3 + 2])
    aa4Coordinates = (coordinates[v3 * 3 + 0], coordinates[v3 * 3 + 1],
                      coordinates[v3 * 3 + 2])
    return (aa1Coordinates, aa2Coordinates, aa3Coordinates, aa4Coordinates)

def getTetraIndicesByAminoAcidIndex(indexList, usedAmino):
    v0 = usedAmino.index(indexList[0])
    v1 = usedAmino.index(indexList[1])
    v2 = usedAmino.index(indexList[2])
    v3 = usedAmino.index(indexList[3])

    return (v0 * 3 + 0, v0 * 3 + 1, v0 * 3 + 2, v1 * 3 + 0, v1 * 3 + 1, v1 * 3 + 2,
            v2 * 3 + 0, v2 * 3 + 1, v2 * 3 + 2, v3 * 3 + 0, v3 * 3 + 1, v3 * 3 + 2,)

# Find the volume of a tetrahedron by the coordinates of its vertices
def findVolumeByCoordinates(((a, b, c), (d, e, f), (g, h, i), (p, q, r))):
    return (1.0 / 6.0) * scipy.linalg.det(scipy.matrix([[a, d, g, p], [b, e, h, q], [c, f, i, r], [1, 1, 1, 1]]))


def calculateJacobian(tetraCoordinates, motifList, usedAmino):

    jMatrix = scipy.zeros((len(motifList), len(tetraCoordinates)))

    currentRowIndex = 0
    for motif in motifList:
        ((x1, x2, x3), (x4, x5, x6),
         (x7, x8, x9), (x10, x11, x12)) = getCoordinatesByAminoAcidIndex(motif["index"],
            tetraCoordinates, usedAmino)
        (dx1loc, dx2loc, dx3loc, dx4loc,
         dx5loc, dx6loc, dx7loc, dx8loc,
         dx9loc, dx10loc, dx11loc, dx12loc) = getTetraIndicesByAminoAcidIndex(motif["index"], usedAmino)

        #import pdb; pdb.set_trace()

        dx1 = (x5 - x11)*(x9 - x12) - (x8 - x11)*(x6 - x12)
        dx2 = (x7 - x10)*(x6 - x12) - (x4 - x10)*(x9 - x12)
        dx3 = (x4 - x10)*(x8 - x11) - (x7 - x10)*(x5 - x11)

        dx4 = (x8 - x11)*(x3 - x12) - (x2 - x11)*(x9 - x12)
        dx5 = (x1 - x10)*(x9 - x12) - (x7 - x10)*(x3 - x12)
        dx6 = (x7 - x10)*(x2 - x11) - (x1 - x10)*(x8 - x11)

        dx7 = (x2 - x11)*(x6 - x12) - (x5 - x11)*(x3 - x12)
        dx8 = (x4 - x10)*(x3 - x12) - (x1 - x10)*(x6 - x12)
        dx9 = (x1 - x10)*(x5 - x11) - (x4 - x10)*(x2 - x11)

        dx10 = -(x8 - x11)*(x3 - x12) + (x5 - x11)*(x3 - x12) \
                - (x5 - x11)*(x9 - x12) + (x2 - x11)*(x9 - x12) \
                + (x8 - x11)*(x6 - x12) - (x2 - x11)*(x6 - x12)
        dx11 = (x4 - x10)*(x9 - x12) - (x1 - x10)*(x9 - x12) \
                - (x7 - x10)*(x6 - x12) + (x1 - x10)*(x6 - x12) \
                + (x7 - x10)*(x3 - x12) - (x4 - x10)*(x3 - x12)
        dx12 = -(x4 - x10)*(x8 - x11) + (x1 - x10)*(x8 - x11) \
                + (x7 - x10)*(x5 - x11) - (x1 - x10)*(x5 - x11) \
                - (x7 - x10)*(x2 - x11) + (x4 - x10)*(x2 - x11)
        currentRow = scipy.zeros(len(tetraCoordinates))

        currentRow[dx1loc] = dx1; currentRow[dx2loc] = dx2; currentRow[dx3loc] = dx3;
        currentRow[dx4loc] = dx4; currentRow[dx5loc] = dx5; currentRow[dx6loc] = dx6;
        currentRow[dx7loc] = dx7; currentRow[dx8loc] = dx8; currentRow[dx9loc] = dx9;
        currentRow[dx10loc] = dx10; currentRow[dx11loc] = dx11; currentRow[dx12loc] = dx12;

        jMatrix[currentRowIndex] = currentRow
        currentRowIndex += 1

    return jMatrix


# Function that is passed into the scipy.optimize function to find the
# possible errors of the tetrahedron volumes
def findResiduals(currentCoordinates, motifList, usedAmino):

    errors = scipy.array([])
    errors = errors.astype(scipy.float64)
    global findResidualsCount

    #yf = tetraCoordinates[0:6]
    #motifCount = 0
    #print "({0:3.3f},{1:3.3f},{2:3.3f}), ({3:3.3f}),{4:3.3f},{5:3.3f})".format(yf[0],yf[1],yf[2], yf[3], yf[4], yf[5])

    for motif in motifList:
        # Extract the volume from the motif.
        measuredVolume = motif["nVol"]
        calculatedVolume = findVolumeByCoordinates(
            getCoordinatesByAminoAcidIndex(motif["index"], currentCoordinates, usedAmino))
        errors = scipy.append(errors, calculatedVolume - measuredVolume)

        #motifCount += 1

        #if motifCount % 10000 == 0: print "mc: ",motifCount

    #if findResidualsCount % 10 == 0:
    #    print "rc", findResidualsCount

    errors = scipy.array(errors, dtype=scipy.float64)

    #print errors
    #raw_input()
    return errors

def findResidualsSquared(tetraCoordinates, motifList, usedAmino):
    errors = findResiduals(tetraCoordinates, motifList, usedAmino)
    return reduce(lambda acc, next: acc + next**2, errors)

#def squareOfResiduals(tetraCoordinates):
#    global listOfMotifsWithVolumes
#    return scipy.sum(scipy.square(findResiduals(tetraCoordinates, listOfMotifsWithVolumes)))

#def cb(yf):
#    print "({0:3.3f},{1:3.3f},{2:3.3f}), ({3:3.3f}),{4:3.3f},{5, 3.3f})".format(yf[0],yf[1],yf[2], yf[3], yf[4], yf[5])


def generateInitialValues(sizeOfProtein):
    x0 = []
    listOfDifferences = \
        generateAminoAcidsDifferenceList()

    currentCoordinates = scipy.array([0,0,0])
    x0.extend(currentCoordinates.tolist())

    for i in range(sizeOfProtein - 1):
        nextIndex = random.randrange(0, len(listOfDifferences))
        diffVector = scipy.array(listOfDifferences[nextIndex])
        currentCoordinates = currentCoordinates + diffVector
        x0.extend(currentCoordinates)

    return x0

def generateInitialValuesFromBackbone(sizeOfProtein, backboneIndices, backboneCoordinates):
    x0 = []
    listOfDifferences = \
        generateAminoAcidsDifferenceList()
    
    # First, instantiate all coordinates to -1 and keep track of
    # how many amino acids have coordinates
    x0 = [-1] * sizeOfProtein
    
    # Iterate through the backbone indices and add the backbone 
    # coordinates
    for i in xrange(len(backboneIndices)):
        x0[backboneIndices[i]-1] = scipy.array(backboneCoordinates[i])
    
    # Iterate through backbone indices again, and this time 
    # attempt to add coordinates to the backbone by traveling forward
    # and traveling backwards  
    
    for index in backboneIndices:
        # Convert index to an array index
        index = index - 1        
        
        # First, iterate backwards
        currentIndex = index - 1
        while currentIndex >= 0 and isinstance(x0[currentIndex], int):
            nextIndex = random.randrange(0, sizeOfProtein)
            diffVector = scipy.array(listOfDifferences[nextIndex])
            x0[currentIndex] = x0[currentIndex + 1] + diffVector
            currentIndex -= 1
        
        # Next, iterate forwards
        currentIndex = index + 1
        while currentIndex < sizeOfProtein and isinstance(x0[currentIndex], int):
            nextIndex = random.randrange(0, sizeOfProtein)
            diffVector = scipy.array(listOfDifferences[nextIndex])
            x0[currentIndex] = x0[currentIndex - 1] + diffVector
            currentIndex += 1
    
    return x0
    

def findCoordinates(motifList, usedAmino):
    x0 = scipy.array(generateInitialValues(len(usedAmino)))

    if Constants.USE_MINIMIZE:
        return scipy.optimize.minimize(findResidualsSquared, x0, args=(motifList, usedAmino),
                                       options={"disp":True})
    elif Constants.USE_OUR_GRADIENT:
        return scipy.optimize.leastsq(findResiduals, x0, args=(motifList, usedAmino),
                                      Dfun=calculateJacobian)[0]
                                      #xtol=0.01, gtol=0.01)
    else:
        return scipy.optimize.leastsq(findResiduals, x0, args=(motifList, usedAmino),
                                      ftol=0.01, xtol=0.01)[0]

def findCoordinatesFromInitial(motifList, usedAmino, initialCoordinatesStraight):
    if Constants.USE_MINIMIZE:
        results =  scipy.optimize.minimize(findResidualsSquared, initialCoordinatesStraight, args=(motifList, usedAmino),
                                       options={"disp":True}, method='Nelder-Mead')
        print results
        return results
    elif Constants.USE_OUR_GRADIENT:
        results = scipy.optimize.leastsq(findResiduals, initialCoordinatesStraight, args=(motifList, usedAmino),
                                Dfun=calculateJacobian)
        print results
        return results[0]
    else:
        results = scipy.optimize.leastsq(findResiduals, initialCoordinatesStraight, args=(motifList, usedAmino),
                                      ftol=0.01, xtol=0.01)
        print results
        return results[0]

def findCoordinatesFromBackbone(motifList, usedAmino, backboneIndices, backboneCoordinates):
    # If len(usedAmino) != sizeOfProtein, this might return an error
    x0 = generateInitialValuesFromBackbone(len(usedAmino), backboneIndices, backboneCoordinates)

    if Constants.USE_MINIMIZE:
        return scipy.optimize.minimize(findResidualsSquared, x0, args=(motifList, usedAmino),
                                       options={"disp":True})
    if Constants.USE_OUR_GRADIENT:
        return scipy.optimize.leastsq(findResiduals, x0, args=(motifList, usedAmino),
                                      Dfun=calculateJacobian)[0]
    else:
        return scipy.optimize.leastsq(findResiduals, x0, args=(motifList, usedAmino), ftol=0.01, xtol=0.01)[0]

def generateCoordinatesFromInitial(motifList, usedAmino, initialCoordinateSet):
    coordinates = []

    for initialCoordinatesGrouped in initialCoordinateSet:
        currentCoordinates = []
        initialCoordinatesStraight = []
        for initialCoordinatesGroup in initialCoordinatesGrouped:
            initialCoordinatesStraight.extend(initialCoordinatesGroup)

        answer = findCoordinatesFromInitial(motifList, usedAmino, initialCoordinatesStraight)
        
        if Constants.USE_MINIMIZE:
            return answer
        
        for i in range(0, len(answer), 3):
            currentCoordinates.append((answer[i], answer[i+1], answer[i+2]))

        coordinates.append(currentCoordinates)

    return coordinates


def generateTenCoordinates(motifList, usedAmino):
    coordinates = []
    for j in range(10):
        currentCoordinates = []
        answer = findCoordinates(motifList, usedAmino)
        for i in range(0, len(answer), 3):
            currentCoordinates.append((answer[i], answer[i+1], answer[i+2]))
        coordinates.append(currentCoordinates)
    return coordinates

def generateTenCoordinatesWithBackbone(motifList, usedAmino, backboneIndices, backboneCoordinates):
    coordinates = []
    for j in range(10):
        currentCoordinates = []
        answer = findCoordinatesFromBackbone(motifList, usedAmino, backboneIndices, backboneCoordinates)
        for i in range(0, len(answer), 3):
            currentCoordinates.append((answer[i], answer[i+1], answer[i+2]))
        coordinates.append(currentCoordinates)
    return coordinates

#def findCoordinates2(sizeOfAminoAcid):
#    x0 = scipy.zeros(sizeOfAminoAcid*3)
#    return scipy.optimize.minimize(squareOfResiduals, x0, callback=cb)
