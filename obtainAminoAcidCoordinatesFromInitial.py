import Constants
import MotifsWithLowestStd
import MotifsLeastSquares
import ReadCoordinatesFile

INITIAL_COORDINATES_FILE = "data/InitialCoordinates-New.txt"
GENERATED_COORDINATES_FILE = "data/GeneratedFromInitial-NewWith1D.txt"

if __name__ == "__main__":
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    initialCoordinateSet = ReadCoordinatesFile.readCoordinateFile(INITIAL_COORDINATES_FILE)["coordSet"]

    # Temporarily only use 30 motifs
    coordinates = MotifsLeastSquares.generateCoordinatesFromInitial(filteredOutput["motifs"][0:200], filteredOutput["visited"],
                                                                    initialCoordinateSet)

    print coordinates
    
    with open(GENERATED_COORDINATES_FILE, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(coordinates) + "\n")
        

    
    
