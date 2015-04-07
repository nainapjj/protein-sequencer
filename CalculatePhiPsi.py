import MMTK.Proteins

PDB_NAME = ""


def calculatePhiPsi(pdbName, output):
    """ Get a list of the phi psi angles of the first chain of the protein """

    protein = MMTK.Proteins.Protein(pdbName, model="no_hydrogens")

    if output:
        for chain in protein:
            print "%s length %i" % (chain.name, len(chain)),
            print "from %s to %s" % (chain[0].name, chain[-1].name)
            for residue in chain :
                print residue.name, residue.phiPsi()
    else:
        for chain in protein:
            # Assume we only have one chain
            phiPsi = []
            for residue in chain:
                phiPsi.append(residue.phiPsi())

            return phiPsi