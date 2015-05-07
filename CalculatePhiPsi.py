import MMTK.Proteins
import Bio.PDB

import re
import glob

import math

PDB_NAME = ""


def calculatePhiPsi(pdbName, output):
    """ Get a list of the phi psi angles of the first chain of the protein """

    configuration = MMTK.PDB.PDBConfiguration(pdbName)
    configuration.deleteHydrogens()
    protein = MMTK.Proteins.Protein(configuration.createPeptideChains(
        model="no_hydrogens"))

    if output:
        for chain in protein:
            print "All angles are in radians"
            print ""
            print "%s length %i" % (chain.name, len(chain)),
            print "from %s to %s" % (chain[0].name, chain[-1].name)
            for residue in chain :
                print residue.name, residue.phiPsi()
    else:
        for chain in protein:
            # Assume we only have one chain
            phiPsi = []
            for residue in chain:
                # Convert the PhiPsi
                phiDeg, psiDeg = None, None
                try:
                    if residue.phiPsi()[0]:
                        phiDeg = math.degrees(residue.phiPsi()[0])
                    if residue.phiPsi()[1]:
                        psiDeg = math.degrees(residue.phiPsi()[1])
                    phiPsi.append((phiDeg, psiDeg))
                except AttributeError:
                    print "Incomplete protein detected. Breaking."
                    return False

            return phiPsi


def findRootMeanSquareDifference(pdbOne, pdbTwo):
    print "Loading PDB file %s" % pdbOne
    structureOne = Bio.PDB.PDBParser().get_structure("abc", pdbOne)
    structureTwo = Bio.PDB.PDBParser().get_structure("abc", pdbTwo)

    print "Everything aligned to first model..."
    ref_model = structureOne[0]
    alt_model = structureTwo[0]

    #Build paired lists of c-alpha atoms, ref_atoms and alt_atoms
    ref_atoms = []
    alt_atoms = []
    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res, index in zip(ref_chain, alt_chain, range(100)):
            assert ref_res.resname == alt_res.resname
            ref_atoms.append(ref_res['CA'])
            alt_atoms.append(alt_res['CA'])

    #Align these paired atom lists:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)

    #Update the structure by moving all the atoms in
    #this model (not just the ones used for the alignment)
    super_imposer.apply(alt_model.get_atoms())

    print "RMS(pdbOne: %s pdbTwo: %s) = %0.2f" % (pdbOne, pdbTwo, super_imposer.rms)

    return super_imposer.rms

# Quick temp script to calculate the root mean square differences
def main_calculateRmsd():
    glob_string = "hamidMethod/*.pdb"
    extract_base_pattern = r"[0-9][A-Z0-9][A-Z0-9][A-Z0-9]"

    for pdbModel in glob.glob(glob_string):
        try:
            baseName = re.search(extract_base_pattern, pdbModel).group(0)
            findRootMeanSquareDifference("pdbs/%s.pdb" % baseName, pdbModel)
        except:
            continue

main_calculateRmsd()
#findRootMeanSquareDifference("pdbs/2RJY.pdb", "hamidMethod/2RJY-gen-tub.rebuilt.pdb")