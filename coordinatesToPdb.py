import ReadCoordinatesFile
import Constants

COORDINATES_FILE = "GeneratedFromInitialNoGradient.txt"
GENERATED_PDB_FILE = "GeneratedFromInitialNoGradient.pdb"

def readAminoAcidFile(aminoFile):
    with open(aminoFile, 'rU') as f:
        return map(str.strip, f.readlines())

def generatePdbFile(coordinateSeq, aminoSeq):
    if len(aminoSeq) < len(coordinateSeq):
        return ()

    pdbString = ""
    currentAmino = 0
    currentCount = 0

    for caCoord in coordinateSeq:
        currentAmino += 1
        #             ATOM   Atom #  Atm Na Res Name ChainId Res #  X Coord Y Coord Z Coord
        pdbString += "ATOM  {0:5d}  {1:2s}  {2:3s} {3:1s}{4:4d}    {5:8.3f}{6:8.3f}{7:8.3f}\n".format(currentAmino, 
                     "CA", aminoSeq[currentCount], "A", currentAmino, caCoord[0], caCoord[1], caCoord[2])
        currentCount += 1

    return pdbString

if __name__ == "__main__":
    firstCoords = ReadCoordinatesFile.readCoordinateFile(COORDINATES_FILE)["coordSet"][0]
    aminSeq = readAminoAcidFile(Constants.aminoFile)
    fileString = generatePdbFile(firstCoords, aminSeq)
    
    with open(GENERATED_PDB_FILE, 'w') as f:
        f.write(fileString)