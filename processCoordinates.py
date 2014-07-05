import ReadCoordinatesFile

COORDINATES_FILE = "RUN_ONCE.txt"

if __name__ == "__main__":
    coordFileOutput = ReadCoordinatesFile.readCoordinateFile(COORDINATES_FILE)    
    print ReadCoordinatesFile.findAverageLengthForEachSet(
        coordFileOutput["coordSet"], coordFileOutput["usedAmino"])