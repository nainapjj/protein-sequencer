# Script is run: python parse_bp.py FILENAME
#   were FILENAME is the file to get the average, std. dev. data from.

# If this script is improted, it provides a function to 
# obtain statistics from a bp volume file

import sys
import numpy as np

def obtainStatsFromFile(fileName):
    # Attempt to go through the file 
    fileData = ""
    with open(fileName) as f:
        fileData = np.array(map(float, f.readlines())) 
    
    return {"mean": np.mean(fileData),
            "std": np.std(fileData)}

if __name__ == "__main__":
    print obtainStatsFromFile(sys.argv[1])