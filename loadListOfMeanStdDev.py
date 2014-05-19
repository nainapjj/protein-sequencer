# Script is run: python loadListOfMeanStdDev.py FILENAME
#   were FILENAME is the file to get the average, std. dev. data from.

# If this script is improted, it provides
# functions to load the list of mean and std. dev
# from Stu's file.

import sys
import re

STD_DEV_PATTERN = "(?<=standard deviation:)[0-9\.]+"
MEAN_PATTERN = "(?<=mean:)[0-9\.]+"
# Matches 7 groups, need 1, 3, 5, and 7  
BASE_PAIR_FINDER = "(?<=\(\')([A-Z])(\', \')([A-Z])(\', \')([A-Z])(\', \')([A-Z])" 

def getStdDevMeanFromFile(filename):
    currentDict = {} 
    with open(filename) as f:
        for line in f.readlines():
            std_match = re.search(STD_DEV_PATTERN, line)
            stdDev = float(std_match.group(0))
            
            mean_match = re.search(MEAN_PATTERN, line)
            mean = float(mean_match.group(0))
            
            baseMotif_match = re.search(BASE_PAIR_FINDER, line)
            baseMotif = baseMotif_match.group(1) + baseMotif_match.group(3) + baseMotif_match.group(5) + baseMotif_match.group(7)
            # Sort the list alphabetically            
            baseMotif = ''.join(sorted(baseMotif))        
            
            currentDict[baseMotif] = (mean, stdDev)
    
    return currentDict
            
if __name__ == "__main__":
    print getStdDevMeanFromFile(sys.argv[1])