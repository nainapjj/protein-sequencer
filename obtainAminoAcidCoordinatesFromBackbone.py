import MotifsWithLowestStd
import Constants
import ReadCoordinatesFile
import MotifsLeastSquares

BACKBONE_COORDINATES_FILE = "500Motifs.txt"
WHICH_BACKBONE = 0
TEN_COORDINATES_FILE_NAME = "CoordinatesFromBackbone.txt"

if __name__ == "__main__":
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    print filteredOutput["visited"]    
    
    backboneFile = ReadCoordinatesFile.readCoordinateFile(BACKBONE_COORDINATES_FILE)

    tenCoordinates = MotifsLeastSquares.generateTenCoordinatesWithBackbone(filteredOutput["motifs"],
                                                                           filteredOutput["visited"], backboneFile["usedAmino"],
                                                                           backboneFile["coordSet"][WHICH_BACKBONE])

    print tenCoordinates
    
    with open(TEN_COORDINATES_FILE_NAME, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(tenCoordinates) + "\n")