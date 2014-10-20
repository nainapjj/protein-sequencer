import scipy
import scipy.linalg
import scipy.optimize
import random
import re
import math

import Constants
import ReadCoordinatesFile

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
    return math.fabs(
        (1.0 / 6.0) * scipy.linalg.det(scipy.matrix([[a, d, g, p], [b, e, h, q], [c, f, i, r], [1, 1, 1, 1]])))


def tub(value, min, max, weight):
    return weight * (scipy.maximum(min - value, 0) + scipy.maximum(value - max, 0))

def groupCoordinatesToXYZTuples(singleCoordinates):
    returnCoordinates = []

    for i in range(len(singleCoordinates) / 3):
        returnCoordinates.append((singleCoordinates[3*i], singleCoordinates[3*i+1],
                                  singleCoordinates[3*i+2]))

    return returnCoordinates

# Function that is passed into the scipy.optimize function to find the
# possible errors of the tetrahedron volumes
def findResiduals(currentCoordinates, motifList, usedAmino):
    errors = scipy.array([])
    errors = errors.astype(scipy.float64)

    for motif in motifList:
        # Extract the volume from the motif.
        measuredVolume = motif["unVol"]
        calculatedVolume = findVolumeByCoordinates(
            getCoordinatesByAminoAcidIndex(motif["index"], currentCoordinates, usedAmino))
        errors = scipy.append(errors, calculatedVolume - measuredVolume)

    if Constants.WHICH_METHOD == Constants.USE_1D_CONSTRAINTS:
        groupedCoordinates = groupCoordinatesToXYZTuples(currentCoordinates)
        lengths = ReadCoordinatesFile.findLengthForCoordinateSet(groupedCoordinates, usedAmino)
        if Constants.USE_TUB:
            lengthResiduals = map(lambda x: tub(x, 0, Constants.MAX_LENGTH, 1000), lengths)
        else:
            lengthResiduals = map(lambda x: (1.0 / Constants.MAX_LENGTH * x) ** 5, lengths)
        errors = scipy.append(errors, lengthResiduals)

    errors = scipy.array(errors, dtype=scipy.float64)

    return errors


def findResidualsSquared(tetraCoordinates, motifList, usedAmino):
    errors = findResiduals(tetraCoordinates, motifList, usedAmino)
    return reduce(lambda acc, next: acc + next ** 2, errors)


def generateInitialValues(sizeOfProtein):
    x0 = []
    listOfDifferences = \
        generateAminoAcidsDifferenceList()

    currentCoordinates = scipy.array([0, 0, 0])
    x0.extend(currentCoordinates.tolist())

    for i in range(sizeOfProtein - 1):
        nextIndex = random.randrange(0, len(listOfDifferences))
        diffVector = scipy.array(listOfDifferences[nextIndex])
        currentCoordinates = currentCoordinates + diffVector
        x0.extend(currentCoordinates)

    return x0


def findCoordinates(motifList, usedAmino):
    print "Finding the initial values..."
    x0 = scipy.array(generateInitialValues(len(usedAmino)))

    print "Performing a least squares analysis..."
    if Constants.WHICH_METHOD == Constants.USE_1D_CONSTRAINTS:
        results = scipy.optimize.leastsq(findResiduals, x0, args=(motifList, usedAmino),
                                         ftol=0.001, xtol=0.001)
        return results[0]

def generateTenCoordinates(motifList, usedAmino):
    coordinates = []
    for j in range(10):
        currentCoordinates = []
        answer = findCoordinates(motifList, usedAmino)
        for i in range(0, len(answer), 3):
            currentCoordinates.append((answer[i], answer[i + 1], answer[i + 2]))
        coordinates.append(currentCoordinates)
    return coordinates
