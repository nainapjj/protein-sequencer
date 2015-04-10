import MMTK.Proteins
import math

math.degrees(.787)

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