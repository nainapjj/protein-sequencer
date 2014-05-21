import Constants
import MotifsWithLowestStd
import MotifsLeastSquares

if __name__ == "__main__":
    listOfMotifsWithVolumes = MotifsWithLowestStd.findListOfMotifsWithVolumes()
    filteredOutput = MotifsWithLowestStd.filterOutUnneededMotifs(
        listOfMotifsWithVolumes, Constants.sizeOfProtein)
    filteredListOfMotifs= filteredOutput["motifs"]
    tenCoordinates = MotifsLeastSquares.generateTenCoordinates(filteredListOfMotifs,
                                            Constants.sizeOfProtein)

    print tenCoordinates
