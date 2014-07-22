import CAscrape
import Constants
import MotifsWithLowestStd
import MotifsLeastSquares

PDB_FILE = "data/3ZOB.pdb"
GENERATED_COORDINATES_FILE = "data/GeneratedFromPdb.txt"

if __name__ == "__main__":
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    initialCoordinateSet = [CAscrape.scrapeCoordinatesFromFile(PDB_FILE)["coordinates"]]

    coordinates = MotifsLeastSquares.generateCoordinatesFromInitial(filteredOutput["motifs"], filteredOutput["visited"],
                                                                    initialCoordinateSet)

    print coordinates

    with open(GENERATED_COORDINATES_FILE, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(coordinates) + "\n")
