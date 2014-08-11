import ReadCoordinatesFile
import coordinatesToPdb
import PdbUtilityFunctions
import glob
import os

coordinatePattern = "outputPdbs/*.txt"
coordinateFiles = glob.glob(coordinatePattern)

for coordinateFile in coordinateFiles:
    coordinateBasename = os.path.basename(coordinateFile)
    coordinateFilename = os.path.splitext(coordinateBasename)[0]
    coordinatesPdbFile = "pdbs/" + coordinateFilename + ".pdb"
    
    print "Generating pdb for coord file: ", coordinateFile    
    
    coordFileOutput = ReadCoordinatesFile.readCoordinateFile(coordinateFile)    
    averageLengthSet = ReadCoordinatesFile.findAverageLengthForEachSet(
        coordFileOutput["coordSet"], coordFileOutput["usedAmino"])
    
    coordsObject = ReadCoordinatesFile.readCoordinateFile(coordinateFile)
    lengths = ReadCoordinatesFile.findAverageLengthForEachSet(coordsObject["coordSet"],
                                                              coordsObject["usedAmino"])
    minLength = 1E10
    minLengthIndex = -1
    for i in xrange(len(lengths)):
        if lengths[i] < minLength:
            minLength = lengths[i]
            minLengthIndex = i
    
    aaSequenceWithIndexDict = PdbUtilityFunctions.getAASequenceWithIndexDictFromPdb(coordinatesPdbFile)
    pdbText = coordinatesToPdb.generatePdbFileFromAminoSeqWithIndex(coordsObject["coordSet"][minLengthIndex],
                                                          coordsObject["usedAmino"],
                                                          aaSequenceWithIndexDict)

    print "Writing to pdb file"    
    
    outputFile = "outputPdbs/" + coordinateFilename + "-alphacarbon-gen.pdb"
    with open(outputFile, 'w') as f:
        f.write(pdbText)
        