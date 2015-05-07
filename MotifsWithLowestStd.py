# Script is run: python obtainLowestVolumes.py

# If this script is imported, it provides
# functions to load the list of mean and std. dev
# from Stu's file.

# Modules that is useful for filtering out unnecessary motifs.

import re
import numpy

import Constants
import StuartFindAllPossibleMotifs
import Indexator

STD_DEV_PATTERN = "(?<=standard deviation:)[0-9\.]+"
MEAN_PATTERN = "(?<=mean:)[0-9\.]+"
# Matches 7 groups, need 1, 3, 5, and 7
BASE_PAIR_FINDER = "(?<=\(\')([A-Z])(\', \')([A-Z])(\', \')([A-Z])(\', \')([A-Z])"


def getStdDevMeanFromFile(filename):
    currentDict = {}
    with open(filename) as f:
        for line in f.readlines():
            std_match = re.search(STD_DEV_PATTERN, line)
            stdDev = float(std_match.group(0))

            mean_match = re.search(MEAN_PATTERN, line)
            mean = float(mean_match.group(0))

            baseMotif_match = re.search(BASE_PAIR_FINDER, line)
            baseMotif = baseMotif_match.group(1) + baseMotif_match.group(3) + baseMotif_match.group(5) + baseMotif_match.group(7)
            # Sort the list alphabetically
            baseMotif = ''.join(sorted(baseMotif))

            currentDict[baseMotif] = {"mean": mean, "std": stdDev}

    return currentDict


def findListOfMotifsWithVolumes(pdbFile):
    motifDict = getStdDevMeanFromFile(Constants.meansAndStdDevFile)
    allMotifs = StuartFindAllPossibleMotifs.findAllPossibleMotifs(pdbFile)

    listOfMotifsWithVolumes = []

    print "Matching each motif with its corresponding un-normalized volume and std. deviation..."
    for motif in allMotifs:
        # Check to see motif is in the dictionary of smallest
        # std. deviations of the file.
        if motif["amino"] in motifDict:
            # Get the unnormalized volumes
            scale = Indexator.indexator(motif["index"])
            unVolume = scale * motifDict[motif["amino"]]["mean"]
            stdDev = motifDict[motif["amino"]]["std"]

            listOfMotifsWithVolumes.append({"index": motif["index"],
                "amino": motif["amino"], "unVol": unVolume, "std": stdDev})
    return listOfMotifsWithVolumes


def filterOutUnneededMotifs(listOfMotifsWithVolumes, sizeOfProtein):
    if Constants.WHICH_CHOOSE_METHOD == Constants.USE_LOWEST_VOL:
        return filterOutUnnededMotifsLowestVol(listOfMotifsWithVolumes, sizeOfProtein)
    elif Constants.WHICH_CHOOSE_METHOD == Constants.USE_LOWEST_STD:
        return filterOutUnneededMotifsLowestStd(listOfMotifsWithVolumes, sizeOfProtein)


def filterOutUnnededMotifsLowestVol(listOfMotifsWithVolumes, sizeOfProtein):

    def filterFunc(motif):
        return numpy.std(motif["index"]) > Constants.STD_THRESHOLD

    print "Sorting by smallest motifs in the protein and then filtering out motifs that are too close..."
    sortedMotifs = sorted(listOfMotifsWithVolumes, key=lambda motif: motif["unVol"])
    filteredMotifs = filter(filterFunc, sortedMotifs)

    return {"motifs": filteredMotifs[:sizeOfProtein * Constants.N_MULTIPLE]}


def filterOutUnneededMotifsLowestStd(listOfMotifsWithVolumes, sizeOfProtein):
    
    newMotifList = []
    unvisitedAA = []
    visitedAA = []

    # First, sort all of the list of motifs by std deviation
    sortedMotifs = sorted(listOfMotifsWithVolumes, key=lambda motif: motif["std"])
    
    # Keep only the first part of the motifs list.
    sortedMotifs = sortedMotifs[:Constants.NUMBER_OF_TOP_MOTIFS]

    # Create a list of size, size of polypeptide
    hasVisitedAA = numpy.zeros(1000) # <-- Temporary fix

    # Create a list to store used motifs and unused motifs
    unusedMotifs = []

    print "Choosing the bottom portion of the motifs to perform least squares on..."
    motifIndex = 0
    print sortedMotifs[0]
    for motif in sortedMotifs:
        try:
            if hasVisitedAA[motif["index"][0]-1] < 9 and hasVisitedAA[motif["index"][1]-1] < 9 and \
               hasVisitedAA[motif["index"][2]-1] < 9 and hasVisitedAA[motif["index"][3]-1] < 9:
                    newMotifList.append(motif)
                    hasVisitedAA[motif["index"][0]-1] += 1
                    hasVisitedAA[motif["index"][1]-1] += 1
                    hasVisitedAA[motif["index"][2]-1] += 1
                    hasVisitedAA[motif["index"][3]-1] += 1
            if(motif["unVol"] < 4.25):
                unusedMotifs.append(motifIndex) 
            else:
                unusedMotifs.append(motifIndex)
        except Exception:
            import pdb; pdb.set_trace()

        motifIndex += 1

    # Pass through the has visited AA list to see which AA have not been visited.
    aaIndex = 1
    for aa in hasVisitedAA:
        if aa == 0:
            unvisitedAA.append(aaIndex)
        else:
            visitedAA.append(aaIndex)
            
        aaIndex += 1

    print "Number of motifs to perform least squares on:", len(newMotifList)

    # Now, pass through the motif list again, add unused motifs so that newMotifList has
    # at least sizeOfProtein*6 (one for each axis * 3) motifs in it
    #for motif in unusedMotifs:
        # First, check if the motifList has sizeOfProtein motifs
    #    if len(newMotifList) >= sizeOfProtein*6:
    #        break
    #    newMotifList.append(sortedMotifs[motif])

    return {"motifs": newMotifList, "unvisited": unvisitedAA, "visited": visitedAA}
