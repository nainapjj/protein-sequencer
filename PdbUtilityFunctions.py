import re


def findLengthOfProteinFromPdb(pdbFile):
    ATOM_PATTERN = r"ATOM"    

    resMax = -1
    resMin = 1E10  
    
    pdbLines = []
    with open(pdbFile) as f:
        pdbLines = f.readlines()
        
    for pdbLine in pdbLines: 
       if re.match(ATOM_PATTERN, pdbLine):
           # Column numbers from 23 to 26 are the 
           # residual sequence number
           currentRes = int(pdbLine[22:26])
           
           if (currentRes > resMax): resMax = currentRes
           if (currentRes < resMin): resMin = currentRes
    
    return resMax - resMin + 1

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
           