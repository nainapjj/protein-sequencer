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

# if __name__ == "__main__":
#     print getStdDevMeanFromFile(sys.argv[1])

def findListOfMotifsWithVolumes():
    motifDict = getStdDevMeanFromFile(Constants.meansAndStdDevFile)
    allMotifs = StuartFindAllPossibleMotifs.findAllPossibleMotifs()

    listOfMotifsWithVolumes = []

    for motif in allMotifs:
        # Check to see motif is in the dictionary of smallest
        # std. deviations of the file.
        if motif["amino"] in motifDict:
            # Get the normalized volumes
            scale = Indexator.indexator(motif["index"])
            nVolume = scale * motifDict[motif["amino"]]["mean"]
            stdDev = motifDict[motif["amino"]]["std"]

            listOfMotifsWithVolumes.append({"index": motif["index"],
                "amino": motif["amino"], "nVol": nVolume, "std": stdDev})
    return listOfMotifsWithVolumes

def filterOutUnneededMotifs(listOfMotifsWithVolumes, sizeOfProtein):
    newMotifList = []
    unvisitedAA = []

    # First, sort all of the list of motifs by std deviation
    sortedMotifs = sorted(listOfMotifsWithVolumes, key=lambda motif: motif["std"])


    # Create a list of size, size of polypeptide
    hasVisitedAA = numpy.zeros(sizeOfProtein)

    # Create a list to store used motifs
    unusedMotifs = []

    motifIndex = 0
    for motif in sortedMotifs:
        if hasVisitedAA[motif["index"][0]-1] < 6 or hasVisitedAA[motif["index"][1]-1] < 6 or \
            hasVisitedAA[motif["index"][2]-1] < 6 or hasVisitedAA[motif["index"][3]-1] < 6:
                newMotifList.append(motif)
                hasVisitedAA[motif["index"][0]-1] += 1
                hasVisitedAA[motif["index"][1]-1] += 1
                hasVisitedAA[motif["index"][2]-1] += 1
                hasVisitedAA[motif["index"][3]-1] += 1
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

    return {"motifs": newMotifList, "unvisited": unvisitedAA}
