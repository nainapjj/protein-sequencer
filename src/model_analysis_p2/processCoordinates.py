from src.model_util_p2 import ReadCoordinatesFile

COORDINATES_FILE = "HamidMethod.txt"

if __name__ == "__main__":
    coordFileOutput = ReadCoordinatesFile.readCoordinateFile(COORDINATES_FILE)

    print "Output average length of each coordinate set:"

    print ReadCoordinatesFile.findAverageLengthForEachSet(
        coordFileOutput["coordSet"], coordFileOutput["usedAmino"])

    for index in range(len(coordFileOutput["coordSet"])):
        print "Output each length for coordinate set #{0}".format(index)
        print ReadCoordinatesFile.findLengthForCoordinateSet(
            coordFileOutput["coordSet"][index], coordFileOutput["usedAmino"])