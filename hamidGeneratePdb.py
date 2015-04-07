import PdbUtilityFunctions
import cPickle
import glob
import os

globString = os.getcwd() + "/pdbs/*.pdb"

PDB_GLOB = glob.glob(globString)
PICKLE_PATH = os.getcwd() + "/hamidMethod/hamidMethod-"
GEN_PATH = os.getcwd() + "/hamidMethod/"


def main():
    with open(PICKLE_FILE, 'r') as f:
        returnObj = cPickle.load(f)

    # Convert the stats file into a normal coordinates file.
    coordinates = []
    for coords in returnObj["stats"]:
        coordinates.append((coords[0]["mean"], coords[1]["mean"], coords[2]["mean"]))

    aminoSeq = PdbUtilityFunctions.findSequenceOfAminoAcids(PDB_FILE)
    pdbString = PdbUtilityFunctions.generatePdbFileFromAminoSeqWithIndex(coordinates, returnObj["visited"], aminoSeq)

    with open(GEN_PDB_FILE, 'w') as f:
        f.write(pdbString)


def tub_make(tub):
    PICKLE_FILE_ret = PICKLE_PATH + tub[61:-4] + "-tub.pickle"
    return PICKLE_FILE_ret

    
def gen_make(tub):
    FILE_ret = GEN_PATH + tub[61:-4] + "-gen-tub.pdb"
    return FILE_ret

    
if __name__ == "__main__":
    for file in PDB_GLOB:
        print "Formatting...", file
        try:
            PDB_FILE = file
            PICKLE_FILE = tub_make(file)
            GEN_PDB_FILE = gen_make(file)        
            main()
        except:
            print file, "failed."
            continue