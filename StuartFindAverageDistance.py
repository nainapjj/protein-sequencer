import cPickle as pickle
import numpy as np

#Unpickle source file
#List = []
#File = open("/Users/haysb/Desktop/PICKLE.txt", "rb")
#filepick = pickle.load(File)
#for i in filepick:
#    List.append(i)

#Define how to find average distance
Templist = []
def LengthSetOld(x):
    TempList = []
    for indice in range(1, (len(x))):
        TempInt = (abs((x[indice][0]) - (x[indice - 1][0])))
        TempInt2 = (abs((x[indice][1]) - (x[indice - 1][1]) - 1))
        TempInt3 = (abs((x[indice][2]) - (x[indice - 1][2]) - 1))
        Vector = np.power(TempInt, 2) + np.power(TempInt2, 2) + np.power(TempInt3, 2)
        Vector = np.power(Vector, 0.5)
        TempList.append(Vector)
    return (np.mean(TempList))

def LengthSet(x, usedAmino):
    TempList = []
    for indice in range(1, (len(usedAmino))):
        if (usedAmino[indice] - usedAmino[indice-1] == 1):
            #print x
            #print (x[usedAmino[indice]][0])
            TempInt = (abs((x[indice][0]) - (x[indice - 1][0])))
            TempInt2 = (abs((x[indice][1]) - (x[indice - 1][1]) - 1))
            TempInt3 = (abs((x[indice][2]) - (x[indice - 1][2]) - 1))
            Vector = np.power(TempInt, 2) + np.power(TempInt2, 2) + np.power(TempInt3, 2)
            Vector = np.power(Vector, 0.5)
            TempList.append(Vector)
    return (np.mean(TempList))

#Find average distance between every set
#RealSet = LengthSet(List)
#Set1 = LengthSet(set1)
#Set2 = LengthSet(set2)
#Set3 = LengthSet(set3)
#Set4 = LengthSet(set4)
#Set5 = LengthSet(set5)
#Set6 = LengthSet(set6)
#Set7 = LengthSet(set7)
#Set8 = LengthSet(set8)
#Set9 = LengthSet(set9)
#Set10 = LengthSet(set10)
#print (Set1, Set2, Set3, Set4, Set5, Set6, Set7, Set8, Set9, Set10, RealSet)
