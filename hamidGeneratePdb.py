import PdbUtilityFunctions
import cPickle

PICKLE_FILE = "hamidMethod/hamidMethod-1DF4-tub.pickle"
PDB_FILE = "pdbs/1DF4.pdb"
GEN_PDB_FILE = "hamidMethod/1DF4-gen-tub.pdb"

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

if __name__ == "__main__":
    main()