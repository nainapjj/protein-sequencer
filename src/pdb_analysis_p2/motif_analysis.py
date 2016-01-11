import Bio.PDB
import Bio.SeqUtils

import numpy as np
import progressbar

import os.path
import glob
import itertools
import io
import pickle
import sqlite3
import math

# -----------------
# This module takes PDB files, extracts the motifs, and calculates the volumes
# of each of those motifs.

# Motifs in this case are defined as combinations of three or four amino acid
# residues from the proteins.
# -----------------

# Constants
ROOT = "./"
PDB_FOLDER = ROOT + "pdbs"
OUTPUT_DB_TETRAHEDRAL = ROOT + "analysis_output/tetrahedral.db"
OUTPUT_DB_TRIANGULAR = ROOT + "analysis_output/triangular.db"

# -----
# SQLLite Type Conversion Functions

def __adapt_list(lst):
    out = io.BytesIO()
    pickle.dump(lst, out)
    out.seek(0)
    return buffer(out.read())


def __convert_list(text):
    out = io.BytesIO(text)
    return pickle.load(out)

# -----

# ----
# Helper Classes

class CarbonAtom():
    def __init__(self, res_name, coord, index):
        self.res = res_name
        self.coord = coord
        self.index = index


class Motif():

    def __init__(self, name, coords, indicies):
        """
        Initialize the parameters inside the Motif class

        :param name: (string)
        :param coords: ([3] numpy array)
        :param indicies: ([3] integers)
        :return:
        """
        self.name = name
        self.coords = coords
        self.indices = indicies
        self.__measure = None
        self.__scaled_measure = None

    @classmethod
    def generate_from_carbon_atoms(cls, carbon_atoms):
        """
        Initialization from carbon atoms

        :param carbon_atoms: [(Bio.PDB.Atom)] list of bio pdb atoms in this motif
        :return:
        """
        # Sort the items based on the residue name
        sorted(carbon_atoms, key=lambda atom: atom.res)

        motif = Motif("", [], [])

        for atom in carbon_atoms:
            motif.name += atom.res
            motif.coords.append(atom.coord)
            motif.indices.append(atom.index)

        return motif

    @classmethod
    def area_from_coords(cls, coords):
        """
        For size-3 motifs, calcuates the area of the 3D triangle formed by the CA atoms.
        Uses Heron's formula.

        Source: http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
        :return: (float) area
        """
        a = np.linalg.norm(coords[0] - coords[1])
        b = np.linalg.norm(coords[0] - coords[2])
        c = np.linalg.norm(coords[1] - coords[2])

        return cls.area_from_lengths((a, b, c))

    @classmethod
    def area_from_lengths(cls, (a, b, c)):
        """
        Uses Heron's formula to calculate the area of the triangle from three sides.

        Source: http://mathworld.wolfram.com/HeronsFormula.html
        :return:
        """
        s = (a + b + c) / 2.0
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        return area

    @classmethod
    def area_scaling_factor(cls, (c1, c2, c3)):
        """
        Calculates the expected area of the triangle given a Markovian polymer model, calculated simply from
        the deltas of the indices.  Assumes a constant bond length B = 1.

        Formula: R_f^2 = (b^2)N. Source: Teroka's Polymer Solutions

        This area can be used to the normalize the motif volumes calculated.
        :return:
        """
        a = math.sqrt(math.fabs(c1 - c2))
        b = math.sqrt(math.fabs(c1 - c3))
        c = math.sqrt(math.fabs(c2 - c3))

        area = cls.area_from_lengths((a, b, c))
        return area

    def calculate_scaled_area(self):
        """
        Scale the area based on the beads-on-a-string model, and return that value
        :return:
        """
        return self.get_measure() / self.area_scaling_factor(self.indices)

    @classmethod
    def volume_from_coords(cls, coords):
        """
        For size-4 motifs, calcuates the area of the 3D tetrahedron formed by the CA atoms.
        :return: (float) volume
        """
        c1 = np.append(coords[0].reshape(3, 1), np.array(1).reshape(1,1), axis=0)
        c2 = np.append(coords[1].reshape(3, 1), np.array(1).reshape(1,1), axis=0)
        c3 = np.append(coords[2].reshape(3, 1), np.array(1).reshape(1,1), axis=0)
        c4 = np.append(coords[3].reshape(3, 1), np.array(1).reshape(1,1), axis=0)

        return 1.0 / 6.0 * abs(np.linalg.det(np.concatenate((c1, c2, c3, c4), axis=1)))

    @classmethod
    def volume_scaling_factor(cls, (c1, c2, c3, c4)):
        """
        Calculates the expected volume of the tetrahedron given a Markovian polymer model, calculated simply from
        the deltas of the indicies.  Assumes a constant bond length B = 1.

        Formula: R_f^2 = (b^2)N. Source: Teroka's Polymer Solutions

        This volume can be used to the normalize the motif volumes calculated.
        :return:
        """
        d12_sq = math.fabs(c1 - c2)
        d13_sq = math.fabs(c1 - c3)
        d14_sq = math.fabs(c1 - c4)
        d23_sq = math.fabs(c2 - c3)
        d24_sq = math.fabs(c2 - c4)
        d34_sq = math.fabs(c3 - c4)

        # Calculate the coordinates of the tetrahedron
        # Source: http://mathforum.org/dr.math/faq/formulas/faq.irreg.tetrahedron.html
        M = np.matrix([[0, d12_sq, d13_sq, d14_sq, 1], [d12_sq, 0, d23_sq, d24_sq, 1], [d13_sq, d23_sq, 0, d34_sq, 1],
                      [d14_sq, d24_sq, d34_sq, 0, 1], [1, 1, 1, 1, 0]])
        volume = (np.linalg.det(M) / 288) ** .5

        assert volume > 0
        return volume

    def calculate_scaled_volume(self):
        """
        Scale the volume based on the beads-on-a-string model, and return that value
        :return:
        """
        return self.get_measure() / self.volume_scaling_factor(self.indices)

    def get_measure(self):
        if self.__measure is None:
            if len(self.indices) == 3:
                self.__measure = self.area_from_coords(self.coords)
            elif len(self.indices) == 4:
                self.__measure = self.volume_from_coords(self.coords)
            else:
                print "Motifs of size %d are not supported" % len(self.indices)
                raise NotImplementedError()

        return self.__measure

    def get_scaled_measure(self):
        if self.__scaled_measure is None:
            if len(self.indices) == 3:
                self.__scaled_measure = self.calculate_scaled_area()
            elif len(self.indices) == 4:
                self.__scaled_measure = self.calculate_scaled_volume()
            else:
                print "Motifs of size %d are not supported" % len(self.indices)
                raise NotImplementedError()

        return self.__scaled_measure
# ---


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
    print("Initalizing the db at (%s)" % db_file)
    conn = sqlite3.connect(db_file, detect_types=sqlite3.PARSE_DECLTYPES)
    sqlite3.register_adapter(list, __adapt_list)
    sqlite3.register_converter("LIST", __convert_list)

    conn.execute("CREATE TABLE data (pdb_name TEXT, motif_name TEXT, measure REAL, scaled_measure REAL, indices LIST, coordinates LIST)")
    conn.commit()

    pdb_parser = Bio.PDB.PDBParser()
    for pdb_file in glob.glob("%s/*.pdb" % pdb_dir):
        print("Begin processing of (%s)..." % pdb_file)
        pdb = pdb_parser.get_structure(pdb_file, pdb_file)

        for polygon in itertools.combinations([residue for residue in pdb.get_residues() if
                                                residue.get_resname() != "HOH"], N):
            atoms = []

            for residue in polygon:
                # Obtain the CA atom
                ca_atoms = [atom for atom in residue.get_atom() if atom.get_id() == "CA"]
                assert len(ca_atoms) == 1
                ca_atom = ca_atoms[0]

                atoms.append(CarbonAtom(res_name=Bio.SeqUtils.seq1(residue.get_resname()), coord=ca_atom.get_coord(),
                                        index=residue.get_id()[1]))

            motif = Motif.generate_from_carbon_atoms(atoms)
            conn.execute("INSERT INTO data (pdb_name, motif_name, measure, scaled_measure, indices, coordinates) VALUES (?, ?, ?, ?, ?, ?)",
                         [pdb_file, motif.name, motif.get_measure(), motif.get_scaled_measure(), motif.indices, motif.coords])
        conn.commit()

        print("End processing of (%s)" % pdb_file)


print ("Begin the processing of triangular motifs...")
analyze_motifs(PDB_FOLDER, OUTPUT_DB_TRIANGULAR, 3)

print("Begin the processsing of tetrahedral motifs...")
analyze_motifs(PDB_FOLDER, OUTPUT_DB_TETRAHEDRAL, 4)







