import obtainAminoAcidCoordinates
import PdbUtilityFunctions
import glob
import os

pdbPattern = "pdbs/*.pdb"
pdbFiles = glob.glob(pdbPattern)

for pdbFile in pdbFiles:
    sizeOfProtein = PdbUtilityFunctions.findLengthOfProteinFromPdb(pdbFile)
    pdbBasename = os.path.basename(pdbFile)
    pdbFilename = os.path.splitext(pdbBasename)[0]
    
    print "Processing pdb file:", pdbFile
    
    # We have to do a bit of monkey patching.
    obtainAminoAcidCoordinates.MotifsLeastSquares.Constants.pdbFile = pdbFile
    obtainAminoAcidCoordinates.MotifsWithLowestStd.Constants.pdbFile = pdbFile
    obtainAminoAcidCoordinates.MotifsWithLowestStd.StuartFindAllPossibleMotifs.Constants.pdbFile = pdbFile
    obtainAminoAcidCoordinates.Constants.pdbFile = pdbFile
    obtainAminoAcidCoordinates.MotifsLeastSquares.Constants.sizeOfProtein = sizeOfProtein
    obtainAminoAcidCoordinates.MotifsWithLowestStd.StuartFindAllPossibleMotifs.Constants.sizeOfProtein = sizeOfProtein
    obtainAminoAcidCoordinates.MotifsWithLowestStd.Constants.sizeOfProtein = sizeOfProtein
    obtainAminoAcidCoordinates.Constants.sizeOfProtein = sizeOfProtein
    
    # Now simply run the main function for obtainAminoAcidCoordinates
    try:
        obtainAminoAcidCoordinates.main("outputPdbs/" + pdbFilename + ".txt")
    except Exception:
        print "Error processing this file. Move on to the next one."

print "All done!"

    


