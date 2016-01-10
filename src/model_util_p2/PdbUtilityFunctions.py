import re

from src.model_generation_p2 import Constants


def returnHydrophobicAminos():
    with open(Constants.hydrophobicFile) as f:
        hAcids = f.readlines()

    return hAcids

def findHydrophobicAcidsNum(pdbFile):
    seq = findSequenceOfAminoAcids(pdbFile)

    with open(Constants.hydrophobicFile) as f:
        hAcids = f.readlines()

    numHydro = 0

    for aa in seq:
        if aa in hAcids: numHydro += 1

    return numHydro


def findSequenceOfAminoAcids(pdbFile):
    ATOM_CA_PATTERN = r"ATOM.+CA"

    with open(pdbFile) as f:
        pdbLines = f.readlines()

    aminoSeq = []
    for pdbLine in pdbLines:
       if re.match(ATOM_CA_PATTERN, pdbLine):
           # Column numbers from 23 to 26 are the
           # residual sequence number
           currentRes = pdbLine[17:20]
           aminoSeq.append(currentRes)

    return aminoSeq

def findLengthOfProteinFromPdb(pdbFile):
    seq = findSequenceOfAminoAcids(pdbFile)
    return len(seq)

def extractCoordinateListFromPdb(pdbFile):
  ATOM_CA_PATTERN = r"ATOM.+CA"  
  coordinates = [] 

  pdbLines = []
  with open(pdbFile) as f:
    pdbLines = f.readlines() 

  for pdbLine in pdbLines: 
       if re.match(ATOM_CA_PATTERN, pdbLine):
           # Column numbers from 31 to 38 are x
           xCoord = float(pdbLine[30:38].strip())
           # Column numbers from 39 to 46 are y
           yCoord = float(pdbLine[38:46].strip())
           # Column numbers from 47 to 54 are z
           zCoord = float(pdbLine[46:54].strip())
           
           coordinates.append((xCoord, yCoord, zCoord))

  return coordinates

def getAASequenceWithIndexDictFromPdb(pdbFile):
    ATOM_PATTERN = r"ATOM"    
    
    aaSequenceWithIndexDict = {}    
    
    pdbLines = []
    with open(pdbFile) as f:
        pdbLines = f.readlines()
    
    lastRes = -1
    for pdbLine in pdbLines: 
       if re.match(ATOM_PATTERN, pdbLine):
           # Column numbers from 23 to 26 are the 
           # residual sequence number
           currentRes = int(pdbLine[22:26])
           
           if lastRes == -1 or currentRes > lastRes:
               currentResName = pdbLine[17:20]
               aaSequenceWithIndexDict[currentRes] = currentResName
               lastRes = currentRes
    
    return aaSequenceWithIndexDict


def generatePdbFileFromAminoSeqWithIndex(coordinateSeq, usedAmino, aminoSeq):

    pdbString = ""
    currentAtom = 0
    currentCount = 0

    for caCoord in coordinateSeq:
        currentAtom += 1

        #             ATOM   Atom #  Atm Na Res Name ChainId Res #  X Coord Y Coord Z Coord
        pdbString += "ATOM  {0:5d}  {1:2s}  {2:3s} {3:1s}{4:4d}    {5:8.3f}{6:8.3f}{7:8.3f}\n".format(currentAtom,
                     "CA", aminoSeq[currentCount], "A", usedAmino[currentCount], caCoord[0], caCoord[1], caCoord[2])
        currentCount += 1

    return pdbString


def getResidualSequence(pdbFile):
    ATOM_PATTERN = r"ATOM"

    with open(pdbFile) as f:
        pdbLines = f.readlines()

    resArray = []

    lastRes = -1
    for pdbLine in pdbLines:
       if re.match(ATOM_PATTERN, pdbLine):
           # Column numbers from 23 to 26 are the
           # residual sequence number
           currentRes = int(pdbLine[22:26])

           if lastRes == -1 or currentRes > lastRes:
               resArray.append(currentRes)
               lastRes = currentRes

    return resArray

def fromResIndToZeroInd(resInd, resSequence):
    return resSequence.index(resInd)