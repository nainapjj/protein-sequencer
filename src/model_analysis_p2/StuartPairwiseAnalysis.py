N_AllResidues = []
N_HydroPhob = []
N_HydroPhil = []

M_AllResidues = []
M_HydroPhob = []
M_HydroPhil = []

import glob
import os
import cPickle as pickle
import numpy as np
import numpy
import itertools
import time
from multiprocessing import Process
from multiprocessing import Pool
import multiprocessing
import random
import resource
import math
import subprocess
import matplotlib.pyplot as plt
import scipy.stats

dict = {'GLY': 'G', 'PRO': 'P', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'AILE': 'I', 'BILE': 'I', 'MET': 'M', 'CYS': 'C', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W', 'HIS': 'H', 'LYS': 'K', 'ARG': 'R', 'GLN': 'Q', 'ASN': 'N', 'GLU': 'E', 'ASP': 'D', 'SER': 'S', 'THR': 'T', 'ATHR': 'T', 'BTHR': 'T', 'ACYS': 'C', 'BCYS': 'C', 'CCYS': 'C', 'ASER': 'S', 'BSER': 'S', 'AVAL': 'V', 'BVAL': 'V', 'AARG': 'R', 'BARG': 'R', 'AASN': 'N', 'BASN': 'N', 
'AGLN': 'Q', 'BGLN': 'Q', 'ALYS': 'K', 'BLYS': 'K'}

hydrodict = {'I': '3.1', 'V': '2.6', 'L': '2.8', 'F': '3.7', 'C': '2.0', 'M': '3.4', 'A': '1.6', 'G': '1.0', 'T': '1.2', 'S': '0.6', 'W': '1.2', 'Y': '-0.7', 'P': '-0.2', 'H': '-3.0', 'E': '-8.2', 'Q': '-4.1', 'D': '-9.2', 'N': '-4.8', 'K': '-8.8', 'R': '-12.3'}

Templist = []
def LengthSet(x):
    TempList = []
    Combo_2_X = itertools.combinations(x, 2)
    for pair in Combo_2_X:
        DifferenceX = abs(pair[1][0] - pair[0][0])
        DifferenceY = abs(pair[1][1] - pair[0][1])
        DifferenceZ = abs(pair[1][2] - pair[0][2])
        Vector = np.power(DifferenceX, 2) + np.power(DifferenceY, 2) + np.power(DifferenceZ, 2)
        Vector = np.power(Vector, 0.5)
        TempList.append(Vector)
    return (TempList)        

def LengthSetAMINO(x, aminolist):
    List = []
    V_List = []
    for item in range(len(x)):
        newitem = (x[item], aminolist[item])
        List.append(newitem)
    Combo_2 = itertools.combinations(List, 2)    
    for pair in Combo_2:
        AminoSum = float(pair[0][1]) + float(pair[1][1])
        if AminoSum > 3:
            pass
        if AminoSum <= 3:
            continue
        DifferenceX = abs(pair[1][0][0] - pair[0][0][0])
        DifferenceY = abs(pair[1][0][1] - pair[0][0][1])
        DifferenceZ = abs(pair[1][0][2] - pair[0][0][2])
        Vector = np.power(DifferenceX, 2) + np.power(DifferenceY, 2) + np.power(DifferenceZ, 2)
        Vector = np.power(Vector, 0.5)
        V_List.append(Vector)     
    return (V_List)   

def LengthSetPAMINO(x, aminolist):
    List = []
    V_List = []
    for item in range(len(x)):
        newitem = (x[item], aminolist[item])
        List.append(newitem)
    Combo_2 = itertools.combinations(List, 2)    
    for pair in Combo_2:
        AminoSum = float(pair[0][1]) + float(pair[1][1])
        if AminoSum > 3:
            continue
        if AminoSum <= 3:
            pass
        DifferenceX = abs(pair[1][0][0] - pair[0][0][0])
        DifferenceY = abs(pair[1][0][1] - pair[0][0][1])
        DifferenceZ = abs(pair[1][0][2] - pair[0][0][2])
        Vector = np.power(DifferenceX, 2) + np.power(DifferenceY, 2) + np.power(DifferenceZ, 2)
        Vector = np.power(Vector, 0.5)
        V_List.append(Vector)     
    return (V_List)   

ProteinName = raw_input("Type protein name.")

#User path information

pdbfiles_dir='hamidMethod/*.pdb'
write_dir='hamidMethod'
no_processors=8

#List and open pdb files, read values
jimbeam =glob.glob(pdbfiles_dir)

arraytotal = []
aminototal = []

for file in jimbeam:
            filename= open(file, "r")
            print file
            array= []
            aminolist = []
            indexlist = []
            STRINGMOD = 0
            while True:
                line=filename.readline()
                if not line:
                    filename.close()
                    break
                line=line.split()
                if (line[0]=='ATOM') :
                    if (line[2]=='CA'):
                        if (line[4]=='A') :
                            indexlist.append(float(line[5]))
                            if line[6].count('') >= 12:
                                if (line[6])[0] is "-":
                                    STRINGMOD = 1
                                line[6] = line[6].replace("-", " ")
                                sequencesplit = str.split(line[6])
                                arraysix = sequencesplit[0] 
                                arrayseven = sequencesplit[1] 
                                try:
                                    arrayeight = sequencesplit[2]
                                except:
                                    arrayeight = line[7]
                                newarrayseven = '-' + arrayseven
                                arrayseven = newarrayseven
                                if STRINGMOD is 1:
                                    newarraysix = '-' + arraysix
                                    arraysix = newarraysix
                                    STRINGMOD = 0
                                sequence = (float(arraysix), float(arrayseven), float(arrayeight))
                                array.append(sequence)
                                    
                            if line[7].count('') >= 12:
                                line[7] = line[7].replace("-", " ")
                                sequencesplit = str.split(line[7])
                                arrayseven = sequencesplit[0]
                                arrayeight = sequencesplit[1]
                                arrayeight = "-" + arrayeight
                                sequence = (float(line[6]), float(arrayseven), float(arrayeight))
                                array.append(sequence)
                            else:
                                sequence= (float(line[6]), float(line[7]), float(line[8]))
                                array.append(sequence)    
                            amino = (line[3])
                            if (amino.count('') is 5):
                                amino= list(amino)
                                amino.remove(amino[0])
                                amino= str(amino[0]) + str(amino[1]) + str(amino[2])
                            aminolist.append(dict[amino])     
                if (line[0]=='MODEL'):
                    if (line[1]=='2'):
                        filename.close()
                        break
            aminototal.append(aminolist)
            arraytotal.append(array)

AllResidues = []
HydroPhob = []
HydroPhil = []

#RealSet = LengthSet(List)
print "Protein", ProteinName
for index in range(2):        
    Array = LengthSet(arraytotal[index])
    HydroAminoList = []
    for subindex in aminototal[index]:
        HydroAminoList.append(hydrodict[subindex])
    ArrayH = LengthSetAMINO(arraytotal[index], HydroAminoList)
    ArrayP = LengthSetPAMINO(arraytotal[index], HydroAminoList)
    AllResidues.append(Array)
    HydroPhob.append(ArrayH)
    HydroPhil.append(ArrayP)
    if index == 0:
        Header = "NATIVE"
    if index == 1:
        Header = "MODEL"
    print Header
    print "Total Size"
    print "Mean: ", numpy.mean(Array), "Median: ", numpy.median(Array)
    print "Hydrophobic Residues"
    print "Mean: ", numpy.mean(ArrayH), "Median: ", numpy.median(ArrayH)
    print "Hydrophylic Residues"
    print "Mean: ", numpy.mean(ArrayP), "Median: ", numpy.median(ArrayP)

x = AllResidues[0]
y = AllResidues[1]

Hx = HydroPhob[0]
Hy = HydroPhob[1]

Px = HydroPhil[0]
Py = HydroPhil[1]

for a in x:
    N_AllResidues.append(a)
    
for a in y:
    M_AllResidues.append(a)

for a in Hx:
    N_HydroPhob.append(a)
    
for a in Hy:
    M_HydroPhob.append(a)

for a in Px:
    N_HydroPhil.append(a)
    
for a in Py:
    M_HydroPhil.append(a)

plt.scatter(x, y)
print scipy.stats.pearsonr(x, y)        
        
