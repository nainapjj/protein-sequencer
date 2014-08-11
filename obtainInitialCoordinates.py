import Constants
import MotifsLeastSquares

TEN_COORDINATES_FILE_NAME = "InitialCoordinates-New.txt"

def groupCoordinatesToXYZTuples(singleCoordinates):
    returnCoordinates = []    
    
    for i in range(len(singleCoordinates) / 3):
        returnCoordinates.append((singleCoordinates[3*i], singleCoordinates[3*i+1], 
                                  singleCoordinates[3*i+2]))
    
    return returnCoordinates

if __name__ == "__main__":                                               
    initialCoordinates = MotifsLeastSquares.generateInitialValues(Constants.sizeOfProtein)
    initialCoordinates = groupCoordinatesToXYZTuples(initialCoordinates)
        
    with open(TEN_COORDINATES_FILE_NAME, 'w') as f:
        f.write(str(range(Constants.sizeOfProtein)) + "\n")
        f.write(str([initialCoordinates]) + "\n")