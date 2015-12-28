import Bio.PDB
import Bio.SeqUtils

import os.path
import glob
import itertools
import io
import pickle
import sqlite3

#-----------------
# This module takes PDB files, extracts the motifs, and calculates the volumes
# of each of those motifs.

# Motifs in this case are defined as combinations of three or four amino acid
# residues from the proteins.
#-----------------

# Constants
OUTPUT_CSV_TETRAHEDRAL = "../analysis_output/tetrahedral.db"
OUTPUT_CSV_TRIANGULAR = "../analysis_output/triangular.db"

#-----
# SQLLite Type Conversion Functions
def __adapt_list(lst):
    out = io.BytesIO()
    pickle.dump(lst, out)
    out.seek(0)
    return buffer(out.read())

def __convert_list(text):
    out = io.BytesIO(text)
    return pickle.load(out)
#-----

#----
# Helper Classes
class CarbonAtom():
    def __init__(self, res_name, coord, index):
        self.res = res_name
        self.coord = coord
        self.index = index

class Motif():
    def __init__(self, carbon_atoms):
        # Sort the items based on the residue name
        sorted(carbon_atoms, key=lambda atom: atom.res)

        self.__measure = None
        self.name = ""
        self.coords = []
        self.indicies = []

        for atom in carbon_atoms:
            self.name += atom.res
            self.coords.append(atom.coord)
            self.indicies.append(atom.index)

    def calculate_area(self):
        """
        For size-3 motifs, calcuates the area of the 3D triangle formed by the CA atoms.
        :return: (float) area
        """
        return 0.0

    def calculate_volume(self):
        """
        For size-4 motifs, calcuates the area of the 3D tetrahedron formed by the CA atoms.
        :return: (float) volume
        """
        return 0.0

    def get_measure(self):
        if self.__measure == None:
            if len(self.indicies) == 3:
                self.__measure = self.calculate_area()
            elif len(self.indicies) == 4:
                self.__measure = self.calculate_volume()
            else:
                print "Motifs of size %d are not supported" % len(self.indicies)
                raise NotImplementedError()

        return self.__measure
#---


def analyze_motifs(pdb_dir, db_file, N):
    """
    Analyzes the pdb files for motifs of length N and calculate their respective
    areas / volumes.  We store that area, index data, and protein name into a SQLLite db.

    :param pdb_dir: (string) the directory containing the pdbs to analyze
    :param db_file: (string) the sqlite3 db filename
    :param N: (int) the size of the motif
    :return: An array of all of the pdbs analyzed successfully
    """

    # Make the output directory path if it doesn't exist
    if os.path.exists(db_file):
        print("Error: The output db (%s) already exists. Please rename the db; otherwise, the database " +
              "contents would be clobbered.") % db_file
        return

    # Initialize the database
    conn = sqlite3.connect(db_file, detect_types=sqlite3.PARSE_DECLTYPES)
    sqlite3.register_adapter(list, __adapt_list)
    sqlite3.register_converter("LIST", __convert_list)

    conn.execute("CREATE TABLE data (motif_name TEXT, measure REAL, indices LIST, coordinates LIST)")
    conn.commit()

    pdb_parser = Bio.PDB.PDBParser()
    for pdb_file in glob.glob("%s/*.pdb" % pdb_dir):
        pdb = pdb_parser.get_structure(pdb_file, pdb_file)

        for triangle in itertools.combinations([residue for residue in pdb.get_residues() if
                                                residue.get_resname() != "HOH"], N):
            atoms = []

            for residue in triangle:
                # Obtain the CA atom
                ca_atoms = [atom for atom in residue.get_atom() if atom.get_id() == "CA"]
                if len(ca_atoms) != 1: import pdb; pdb.set_trace()
                assert len(ca_atoms) == 1
                ca_atom = ca_atoms[0]

                atoms.append(CarbonAtom(res_name=Bio.SeqUtils.seq1(residue.get_resname()), coord=ca_atom.get_coord(),
                                        index=residue.get_id()[1]))

            motif = Motif(atoms)
            conn.execute("INSERT INTO data (motif_name, measure, indices, coordinates) VALUES (?, ?, ?, ?)",
                         [motif.name, motif.get_measure(), motif.indicies, motif.coords])
            conn.commit()

analyze_motifs("../pdb_temp", "test.db", 3)








