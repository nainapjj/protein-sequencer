import constants
import re

# Groups: 1 - X, 3 - Y, 5 - Z
DIFFERENCE_COORDINATES_PTN = r"(?<=\()([\-\.0-9]+)(, )([\-\.0-9]+)(, )([\-\.0-9]+)"

def generateAminoAcidsDifferenceList():
    listOfAminoDifferences = []
    with open(constants.differenceInAminoFile) as f:
        for line in f.readlines():
            coord_match = re.search(DIFFERENCE_COORDINATES_PTN, line)
            x = float(coord_match.group(1))
            y = float(coord_match.group(3))
            z = float(coord_match.group(5))
            listOfAminoDifferences.append((x, y, z))
    return listOfAminoDifferences