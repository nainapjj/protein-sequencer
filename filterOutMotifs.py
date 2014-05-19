# Modules that is useful for filtering out unnecessary motifs.

import numpy

def filterOutUnneededMotifs(listOfMotifsWithVolumes, sizeOfProtein):
    newMotifList = []    
    unvisitedAA = []
    
    # First, sort all of the list of motifs by std deviation
    sortedMotifs = sorted(listOfMotifsWithVolumes, key=lambda motif: motif[3])
    
    
    
    # Create a list of size, size of polypeptide
    hasVisitedAA = numpy.zeros(sizeOfProtein)
    
    # Create a list to store used motifs
    unusedMotifs = []
    
    motifIndex = 0
    for motif in sortedMotifs:
        if hasVisitedAA[motif[0][0]-1] < 6 or hasVisitedAA[motif[0][1]-1] < 6 or \
            hasVisitedAA[motif[0][2]-1] < 6 or hasVisitedAA[motif[0][3]-1] < 6:
                newMotifList.append(motif)
                hasVisitedAA[motif[0][0]-1] += 1
                hasVisitedAA[motif[0][1]-1] += 1
                hasVisitedAA[motif[0][2]-1] += 1
                hasVisitedAA[motif[0][3]-1] += 1
        else:
            unusedMotifs.append(motifIndex)
            
        motifIndex += 1
    
    # Pass through the has visited AA list to see which AA have not been visited.
    aaIndex = 0
    for aa in hasVisitedAA:
        if aa == 0: unvisitedAA.append(aaIndex)
        aaIndex += 1
    
    # Now, pass through the motif list again, add unused motifs so that newMotifList has
    # at least sizeOfProtein*3 (one for each axis) motifs in it
    for motif in unusedMotifs:
        # First, check if the motifList has sizeOfProtein motifs
        if len(newMotifList) >= sizeOfProtein*6:
            break
        newMotifList.append(sortedMotifs[motif])
    
    return (newMotifList, unvisitedAA)
    