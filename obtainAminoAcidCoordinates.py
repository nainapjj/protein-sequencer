import Constants
import MotifsWithLowestStd
import MotifsLeastSquares

TEN_COORDINATES_FILE_NAME = "500Motifs.txt"

def main(tenCoordinatesFileName):
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    print "Visited AA indices: ", filteredOutput["visited"]    
    tenCoordinates = MotifsLeastSquares.generateTenCoordinates(filteredOutput["motifs"],
                                            filteredOutput["visited"])

    #print tenCoordinates
    
    with open(tenCoordinatesFileName, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(tenCoordinates) + "\n")

if __name__ == "__main__":
    main(TEN_COORDINATES_FILE_NAME)
    
    
