import os.path
import re

import Constants
import ExecuteCommands
from ..model_analysis_p2 import CalculatePhiPsi
from ..model_util_p2 import PdbUtilityFunctions

TEMP_XYZ_OUTPUT_NAME = "outputFile"
TINKER_OUTPUT_PATTERN = r"(.*)_[0-9]+"

PULCHRA_OUTPUT_PATTERN = r"\.rebuilt\.pdb"

FROM_PDB = "hamidMethod/1DF4-gen.pdb"
TO_PDB = "1DF4-gen-tinker 3.pdb"
TO_PDB_2 = "1DF4-gen-tinker 4.pdb"


def __pulchraOutputFileName(originalFile):
    baseName = os.path.splitext(originalFile)[0]
    ext = os.path.splitext(originalFile)[1]

    return "%s.rebuilt%s" % (baseName, ext)


def pulchraReconstructProtein(pdbFile):
    pdbOut = __pulchraOutputFileName(pdbFile)

    # Don't Pulchra pdbs with PDBs that are already generated
    # to prevent possible data loss and save time
    if os.path.exists(pdbOut):
        print "%s already exists, no Pulchra ran" % pdbOut
        return pdbOut

    ExecuteCommands.runCommandWithInput("%s %s" % (Constants.PULCHRA_EXEC, pdbFile), "")

    if os.path.exists(pdbOut):
        return pdbOut
    else:
        print "Pulchra error with %s" % pdbFile
        return False


def __tinkerOutputFileName(originalOutput):
    found = re.search(TINKER_OUTPUT_PATTERN, originalOutput)

    if found:
        baseName = found.group(1)
    else:
        baseName = originalOutput

    number = 2
    newFile = "%s_%d" % (baseName, number)

    while os.path.exists(newFile):
        number += 1
        newFile = "%s_%d" % (baseName, number)

    return newFile


def __tinkerNewOutputFileName(originalOutput, ext):
    baseName = os.path.basename(originalOutput)
    fName = os.path.splitext(baseName)[0]
    newFile = "%s.%s" % (fName, ext)
    number = 2

    while os.path.exists(newFile):
        newFile = "%s.%s_%d" % (fName, ext, number)
        number += 1

    return newFile


def tinkerCreateProtein(xyzOutputFile, aaList, phiPsiAngles):
    # Figure out what the newFile name will be
    xyzFile = __tinkerNewOutputFileName(xyzOutputFile, "xyz")
    intFile = __tinkerNewOutputFileName(xyzOutputFile, "int")
    seqFile = __tinkerNewOutputFileName(xyzOutputFile, "seq")

    # First, construct the input for the Tinker

    # The Output file name, title of the XYZ File, Potential Parameter Filename
    inputString = "%s\n%s\n%s\n" % (xyzOutputFile, xyzOutputFile,
                                    Constants.TINKER_POTENTIAL_PARAMETER)

    # Each of the Amino Acid lists, plus their acid potentials
    for index in xrange(len(aaList)):
        phi, psi = 0, 0
        if phiPsiAngles[index][0]:
            phi = phiPsiAngles[index][0]
        if phiPsiAngles[index][1]:
            psi = phiPsiAngles[index][1]

        inputString += "%s %d %d\n" % (aaList[index], phi, psi)

    # Input ends with two newlines
    inputString += "\n\n"

    ExecuteCommands.runCommandWithInput(Constants.TINKER_PROTEIN_EXEC, inputString)
    return xyzFile, intFile, seqFile

def runEnergyMinimization(xyzInputFile):
    # Figure out new file name
    newFile = __tinkerOutputFileName(xyzInputFile)

    # Construct the input for Tinker

    # First is the xyz input file, potential parameter filename, then default RMS gradient
    inputString = "%s\n%s\n\n" % (xyzInputFile, Constants.TINKER_POTENTIAL_PARAMETER)

    ExecuteCommands.runCommandWithInput(Constants.TINKER_ENERGY_MINIMIZATION, inputString)
    return newFile


def runXyzToPdb(xyzInputFile):
    # Figure out new filename
    newFile = __tinkerNewOutputFileName(xyzInputFile, "pdb")

    # Construct the input for Tinker

    # First is the xyz filename, then potential parameter filename
    inputString = "%s\n%s\n" % (xyzInputFile, Constants.TINKER_POTENTIAL_PARAMETER)

    ExecuteCommands.runCommandWithInput(Constants.TINKER_TO_PDB, inputString)
    return newFile


def main(pdbInputFile, pdbOutputFileMini, pdbOutputFileSame):
    baseName = os.path.splitext(os.path.basename(pdbInputFile))[0]

    # Don't rebuild pdbs that generated from Pulchra
    if re.search(PULCHRA_OUTPUT_PATTERN, pdbInputFile):
        print "Won't run %s directly, is a Pulchra file" % pdbInputFile
        return

    pulchraPdb = pulchraReconstructProtein(pdbInputFile)
    if not pulchraPdb:
        pulchraPdb = pdbInputFile

    aaSeq = PdbUtilityFunctions.findSequenceOfAminoAcids(pulchraPdb)
    phiPsi = CalculatePhiPsi.calculatePhiPsi(pulchraPdb, False)

    if not phiPsi:
        return

    xyzOutput_1, intOutput_1, seqOutput_1 = tinkerCreateProtein(baseName + TEMP_XYZ_OUTPUT_NAME, aaSeq, phiPsi)
    pdbOutputSame = runXyzToPdb(xyzOutput_1)
    xyzOutput_2 = runEnergyMinimization(xyzOutput_1)
    pdbOutputMini = runXyzToPdb(xyzOutput_2)

    # Rename the PDB to desired output
    os.rename(pdbOutputMini, pdbOutputFileMini)
    os.rename(pdbOutputSame, pdbOutputFileSame)

    # Clean up the temp files
    os.remove(xyzOutput_1)
    os.remove(intOutput_1)
    os.remove(seqOutput_1)
    os.remove(xyzOutput_2)


if __name__ == "__main__":
    main(FROM_PDB, TO_PDB, TO_PDB_2)