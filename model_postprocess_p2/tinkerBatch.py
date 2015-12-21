import glob
import multiprocessing
import os.path

from model_generation_p2 import Constants
from model_postprocess_p2 import tinkerReconstructAndEnergyMinimize

GLOB_STRING = "hamidMethod/*.pdb"
TINKER_DIRECTORY = "tinker/"


def runTinker(pdbFile):
    baseFile = os.path.splitext(os.path.basename(pdbFile))[0]
    outputPdbMini = TINKER_DIRECTORY + "%s-tinker-%s.pdb" % (baseFile, Constants.TINKER_ENERGY_MINIMIZATION)
    outputPdbSame = TINKER_DIRECTORY + "%s-tinker-same.pdb" % baseFile

    print "Minimizing %s to %s" % (pdbFile, outputPdbMini)
    print "Translating %s to %s" % (pdbFile, outputPdbSame)

    tinkerReconstructAndEnergyMinimize.main(pdbFile, outputPdbMini, outputPdbSame)


def main():
    pdbGlob = glob.glob(GLOB_STRING)

    if pdbGlob:
        print "You have %d virtual CPUs" % multiprocessing.cpu_count()
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(runTinker, pdbGlob)

if __name__ == "__main__":
    main()