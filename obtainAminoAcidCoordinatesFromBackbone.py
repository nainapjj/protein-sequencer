# -*- coding: utf-8 -*-

if __name__ == "__main__":
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    print filteredOutput["visited"]    
    tenCoordinates = MotifsLeastSquares.generateTenCoordinatesWithBackbone(filteredOutput["motifs"],
                                            filteredOutput["visited"])

    print tenCoordinates
    
    with open(TEN_COORDINATES_FILE_NAME, 'w') as f:
        f.write(str(filteredOutput["visited"]) + "\n")
        f.write(str(tenCoordinates) + "\n")