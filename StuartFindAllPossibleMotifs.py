import numpy as np
import numpy
import itertools

import Constants

dict = {'GLY': 'G', 'PRO': 'P', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'AILE': 'I', 'BILE': 'I', 'MET': 'M', 'CYS': 'C', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W', 'HIS': 'H', 'LYS': 'K', 'ARG': 'R', 'GLN': 'Q', 'ASN': 'N', 'GLU': 'E', 'ASP': 'D', 'SER': 'S', 'THR': 'T', 'ATHR': 'T', 'BTHR': 'T', 'ACYS': 'C', 'BCYS': 'C', 'CCYS': 'C', 'ASER': 'S', 'BSER': 'S', 'AVAL': 'V', 'BVAL': 'V', 'AARG': 'R', 'BARG': 'R', 'AASN': 'N', 'BASN': 'N',
'AGLN': 'Q', 'BGLN': 'Q', 'ALYS': 'K', 'BLYS': 'K'}

pitchdict = {'G': 1, 'P': 2, 'A': 3, 'V': 4, 'L': 5, 'I': 6, 'M': 7, 'C': 8, 'F': 9, 'Y': 10, 'W': 11, 'H': 12, 'K': 13, 'R': 14, 'Q': 15, 'N': 16, 'E': 17, 'D': 18, 'S': 19, 'T': 20}

def indexator((a, b, c, d)):
    u = abs(a - d)
    v = abs(a - c)
    w = abs(a - b)
    U = abs(b - c)
    V = abs(b - d)
    W = abs(c - d)
    M= numpy.matrix([[0, u, v, w, 1], [u, 0, W, V, 1], [v, W, 0, U, 1], [w, V, U, 0, 1], [1, 1, 1, 1, 0]])
    volume = (numpy.linalg.det(M)/288)**.5
    M=[]
    volume = round(volume, 7)
    return volume

def ParVol((a,  b,  c,  d)):
    listing = [a,  b,  c,  d]
    listing = list(listing)
    o = min(listing)
    listing.remove(o)
    point1 = listing[0] - o
    point2 = listing[1] - o
    point3 = listing[2] - o
    newnum = ((1.0/3.0) * (np.power((point1), 0.5) * np.power((point2),  0.5) * np.power((point3),  0.5)))
    return newnum

# Creates a list with FORMAT: ((ind1, ind2, ind3, ind4), "am1am2am3am4")
def findAllPossibleMotifs():
    candidate = Constants.pdbFile
    filename= open(candidate, "r")
    datalist = []

    # Reads the PDB file.
    while True:
        line=filename.readline()
        if not line:
            filename.close()
            break
        line=line.split()
        if (line[0]=='ATOM') :
            if (line[2]=='CA'):
                if (line[4]=='A'):
                    amino = dict[line[3]]
                    index = float(line[5])
                    #amino = pitchdict[amino]
                    datalist.append((amino, index))
        if line[0]=='ENDMDL':
            break
    filename.close()

    # Generates all of the combinations of Amino Acids (to find all
    # possible motifs)
    datacombo = list(itertools.combinations(datalist, 4))

    # Creates a list with FORMAT: ({"index": ind1, ind2, ind3, ind4),
    # "amino": "am1am2am3am4"})
    dataformat = []
    for i in datacombo:
        index = i[0][1], i[1][1], i[2][1], i[3][1]
        amino = i[0][0] + i[1][0] + i[2][0] + i[3][0]
        amino = ''.join(sorted(amino))
        dataformat.append({"index": index, "amino": amino})

    return dataformat
