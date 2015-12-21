import glob
import numpy as np
import numpy
import itertools

dict = {'GLY': 'G', 'PRO': 'P', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'AILE': 'I', 'BILE': 'I', 'MET': 'M', 'CYS': 'C', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W', 'HIS': 'H', 'LYS': 'K', 'ARG': 'R', 'GLN': 'Q', 'ASN': 'N', 'GLU': 'E', 'ASP': 'D', 'SER': 'S', 'THR': 'T', 'ATHR': 'T', 'BTHR': 'T', 'ACYS': 'C', 'BCYS': 'C', 'CCYS': 'C', 'ASER': 'S', 'BSER': 'S', 'AVAL': 'V', 'BVAL': 'V', 'AARG': 'R', 'BARG': 'R', 'AASN': 'N', 'BASN': 'N', 
'AGLN': 'Q', 'BGLN': 'Q', 'ALYS': 'K', 'BLYS': 'K'}

#Function takes 4 x, y, z coordinates, finds volume in the form of a tetrahedron.        
def side(((xa, ya, za), (xb, yb, zb), (xc, yc, zc), (xd, yd, zd))):
    vectorAa= (xa-xb)
    vectorAb= (ya-yb)
    vectorAc= (za-zb)
    vectorBa= (xa-xc)
    vectorBb= (ya-yc)
    vectorBc= (za-zc)
    vectorCa= (xb-xc)
    vectorCb= (yb-yc)
    vectorCc= (zb-zc)
    vectorDa= (xa-xd)
    vectorDb= (ya-yd)
    vectorDc= (za-zd)
    vectorEa= (xb-xd)
    vectorEb= (yb-yd)
    vectorEc= (zb-zd)
    vectorFa= (xc-xd)
    vectorFb= (yc-yd)
    vectorFc= (zc-zd)
    baseA= float(np.power(vectorAa, 2) + np.power(vectorAb, 2))
    U= float(baseA + np.power(vectorAc, 2))
    baseB= float(np.power(vectorBa, 2) + np.power(vectorBb, 2))
    v= float(baseB + np.power(vectorBc, 2))
    baseC= float(np.power(vectorCa, 2) + np.power(vectorCb, 2))
    w= float(baseC + np.power(vectorCc, 2))
    baseD= float(np.power(vectorDa, 2) + np.power(vectorDb, 2))
    W= float(baseD + np.power(vectorDc, 2))
    baseE= float(np.power(vectorEa, 2) + np.power(vectorEb, 2))
    V= float(baseE + np.power(vectorEc, 2))
    baseF= float(np.power(vectorFa, 2) + np.power(vectorFb, 2))
    u= float(baseF + np.power(vectorFc, 2))
    M= numpy.matrix([[0, u, v, w, 1], [u, 0, W, V, 1], [v, W, 0, U, 1], [w, V, U, 0, 1], [1, 1, 1, 1, 0]])
    #Z= M
    volume = (numpy.linalg.det(M)/288)**.5
    #Z=[]
    M=[]
    return volume


def scrapeCoordinatesFromFile(pdbfiles_dir):
    #List and open pdb files, read values
    jimbeam =glob.glob(pdbfiles_dir)
    for file in jimbeam:
        filename= open(file, "r")
        array= []
        aminolist = []
        full = []
        STRINGMOD = 0
        while True:
            line=filename.readline()
            if not line:
                filename.close()
                break
            line=line.split()
            if (line[0]=='HEADER'):
                directoryname= line[-1]
                print "Reading",  directoryname
            if (line[0]=='ATOM') :
                if (line[2]=='CA'):
                    if (line[4]=='A') or (line[4]=='X'):
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
                        full.append(line[0:6])
            if (line[0]=='MODEL'):
                if (line[1]=='2'):
                    filename.close()
                    break

    return { "aminoName": aminolist, "coordinates": array }