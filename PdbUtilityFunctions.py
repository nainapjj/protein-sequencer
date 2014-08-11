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
           