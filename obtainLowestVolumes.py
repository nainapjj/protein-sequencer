# Script is run: python obtainLowestVolumes.py

import constants
import loadListOfMeanStdDev
import findAllPossibleMotifs
import indexator

if __name__ == "__main__":
    motifDict = loadListOfMeanStdDev.getStdDevMeanFromFile(constants.meansAndStdDevFile)
    allMotifs = findAllPossibleMotifs.findAllPossibleMotifs()
    
    listOfMotifsWithVolumes = []    
    
    for motif in allMotifs:
        # Check to see motif is in the dictionary of smallest
        # std. deviations of the file.
        if motif[1] in motifDict:
            # Get the normalized volumes
            scale = indexator.indexator(motif[0])
            nVolume = scale * motifDict[motif[1]][0]
            stdDev = motifDict[motif[1]][1]
            
            listOfMotifsWithVolumes.append((motif[0], motif[1], nVolume, stdDev))
            
