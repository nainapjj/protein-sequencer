import MotifsWithLowestStd
import MotifsLeastSquares
import PdbUtilityFunctions

PDB_FILE = "data/3ZOB-one.pdb"
TEN_COORDINATES_FILE_NAME = "HamidMethod.txt"

def main(tenCoordinatesFileName):
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes(PDB_FILE)
    sizeOfProtein = PdbUtilityFunctions.findLengthOfProteinFromPdb(PDB_FILE)

    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, sizeOfProtein)
    print "Visited AA indices: ", filteredOutput["visited"]    
    tenCoordinates = MotifsLeastSquares.generateTenCoordinates(filteredOutput["motifs"],
                                            filteredOutput["visited"])


    with open(tenCoordinatesFileName, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(tenCoordinates) + "\n")

if __name__ == "__main__":
    main(TEN_COORDINATES_FILE_NAME)
    
    
