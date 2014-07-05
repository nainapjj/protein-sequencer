import Constants
import MotifsWithLowestStd
import MotifsLeastSquares
import ReadCoordinatesFile

INITIAL_COORDINATES_FILE = "InitialCoordinates.txt"
GENERATED_COORDINATES_FILE = "GeneratedFromInitialFromMinimum.txt"

if __name__ == "__main__":
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    initialCoordinateSet = ReadCoordinatesFile.readCoordinateFile(INITIAL_COORDINATES_FILE)["coordSet"]

    coordinates = MotifsLeastSquares.generateCoordinatesFromInitial(filteredOutput["motifs"], filteredOutput["visited"],
                                                                    initialCoordinateSet)

    print coordinates
    
    with open(GENERATED_COORDINATES_FILE, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(coordinates) + "\n")
    
    
