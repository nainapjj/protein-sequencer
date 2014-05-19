import scipy
import scipy.linalg
import scipy.optimize
import random

import loadListOfAminoDifferences

findResidualsCount = 0

# Attempt to obtain the least squares function 
# the matrix in question.
def getCoordinatesByAminoAcidIndex(indexList, tetraCoordinates):
    v0 = int(indexList[0] - 1)
    v1 = int(indexList[1] - 1)
    v2 = int(indexList[2] - 1)
    v3 = int(indexList[3] - 1)
    aa1Coordinates = (tetraCoordinates[v0 * 3 + 0], tetraCoordinates[v0 * 3 + 1], 
                      tetraCoordinates[v0 * 3 + 2])
    aa2Coordinates = (tetraCoordinates[v1 * 3 + 0], tetraCoordinates[v1 * 3 + 1], 
                      tetraCoordinates[v1 * 3 + 2])
    aa3Coordinates = (tetraCoordinates[v2 * 3 + 0], tetraCoordinates[v2 * 3 + 1], 
                      tetraCoordinates[v2 * 3 + 2])
    aa4Coordinates = (tetraCoordinates[v3 * 3 + 0], tetraCoordinates[v3 * 3 + 1], 
                      tetraCoordinates[v3 * 3 + 2])
    return (aa1Coordinates, aa2Coordinates, aa3Coordinates, aa4Coordinates)

def getTetraIndicesByAminoAcidIndex(indexList):
    v0 = indexList[0] - 1
    v1 = indexList[1] - 1
    v2 = indexList[2] - 1
    v3 = indexList[3] - 1
    
    return (v0 * 3 + 0, v0 * 3 + 1, v0 * 3 + 2, v1 * 3 + 0, v1 * 3 + 1, v1 * 3 + 2, 
            v2 * 3 + 0, v2 * 3 + 1, v2 * 3 + 2, v3 * 3 + 0, v3 * 3 + 1, v3 * 3 + 2,)

# Find the volume of a tetrahedron by the coordinates of its vertices
def findVolumeByCoordinates(((a, b, c), (d, e, f), (g, h, i), (p, q, r))):
    return (1.0 / 6.0) * scipy.linalg.det(scipy.matrix([[a, d, g, p], [b, e, h, q], [c, f, i, r], [1, 1, 1, 1]]))
    

def calculateJacobian(tetraCoordinates, motifList):
    
    jMatrix = scipy.zeros((len(motifList), len(tetraCoordinates))) 
    
    currentRowIndex = 0
    for motif in motifList:
        ((x1, x2, x3), (x4, x5, x6), 
         (x7, x8, x9), (x10, x11, x12)) = getCoordinatesByAminoAcidIndex(motif[0],
            tetraCoordinates)
        (dx1loc, dx2loc, dx3loc, dx4loc,
         dx5loc, dx6loc, dx7loc, dx8loc,
         dx9loc, dx10loc, dx11loc, dx12loc) = getTetraIndicesByAminoAcidIndex(motif[0])
        
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
def findResiduals(tetraCoordinates, motifList):  
    
    #print tetraCoordinates
    
    errors = scipy.array([])
    errors = errors.astype(scipy.float64)
    global findResidualsCount
    
    yf = tetraCoordinates[0:6] 
    #motifCount = 0
    #print "({0:3.3f},{1:3.3f},{2:3.3f}), ({3:3.3f}),{4:3.3f},{5:3.3f})".format(yf[0],yf[1],yf[2], yf[3], yf[4], yf[5])
    
    for motif in motifList:
        # Extract the volume from the motif.
        measuredVolume = motif[2]
        calculatedVolume = findVolumeByCoordinates(
            getCoordinatesByAminoAcidIndex(motif[0], tetraCoordinates))
        errors = scipy.append(errors, calculatedVolume - measuredVolume)
        
        #motifCount += 1
        
        #if motifCount % 10000 == 0: print "mc: ",motifCount

    #if findResidualsCount % 10 == 0:
    #    print "rc", findResidualsCount
    
    errors = scipy.array(errors, dtype=scipy.float64)
    
    #print errors
    #raw_input()
    
    return errors

#def squareOfResiduals(tetraCoordinates):
#    global listOfMotifsWithVolumes
#    return scipy.sum(scipy.square(findResiduals(tetraCoordinates, listOfMotifsWithVolumes)))

#def cb(yf):
#    print "({0:3.3f},{1:3.3f},{2:3.3f}), ({3:3.3f}),{4:3.3f},{5, 3.3f})".format(yf[0],yf[1],yf[2], yf[3], yf[4], yf[5])


def generateInitialValues(sizeOfProtein):
    x0 = []
    listOfDifferences = \
        loadListOfAminoDifferences.generateAminoAcidsDifferenceList()
    
    currentCoordinates = scipy.array([0,0,0])
    x0.extend(currentCoordinates.tolist())
    
    for i in range(sizeOfProtein - 1):
        nextIndex = random.randrange(0, sizeOfProtein)
        diffVector = scipy.array(listOfDifferences[nextIndex])
        currentCoordinates = currentCoordinates + diffVector
        x0.extend(currentCoordinates)
    
    return x0

def findCoordinates(motifList, sizeOfAminoAcid):
    #scipy.set_printoptions(threshold=scipy.nan)
    #x0 = ([0.0] * 3)
    #x0.extend([1.0] * ((sizeOfAminoAcid*3) - 3))

    # Generate the random initial value conditions for the 
    # motifs.
    #x0 = []
    #for i in range(sizeOfAminoAcid):
    #    x0.extend([i]*3)
    #x0 = scipy.arange(sizeOfAminoAcid*3)**2
    
    x0 = scipy.array(generateInitialValues(sizeOfAminoAcid))
    
    return scipy.optimize.leastsq(findResiduals, x0, args=(motifList),
                                  Dfun=calculateJacobian)[0]
                                  #xtol=0.01, gtol=0.01)
                                  
def generateTenCoordinates(motifList, sizeOfAminoAcid):
    coordinates = []    
    for j in range(10):
        currentCoordinates = []
        answer = findCoordinates(motifList, sizeOfAminoAcid)
        for i in range(0, len(answer), 3):
            currentCoordinates.append((answer[i], answer[i+1], answer[i+2]))
        coordinates.append(currentCoordinates)
    return coordinates


#def findCoordinates2(sizeOfAminoAcid):
#    x0 = scipy.zeros(sizeOfAminoAcid*3)
#    return scipy.optimize.minimize(squareOfResiduals, x0, callback=cb)