import PdbUtilityFunctions
import itertools
import MotifsLeastSquares
import scipy as sp
import multiprocessing
import math

PDB_FILE_NATIVE = "data/3ZOB-one.pdb"
PDB_FILE_GENERATED = "HamidMethod.pdb"

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def aminoAcidCoordinateMotifs(coordList):
	coordinateGenerator = itertools.combinations(coordList, 4)
	for coordinate in coordinateGenerator:
		yield coordinate

def volumeOfMotifs(coordList):
	numberCombos = nCr(len(coordList), 4)
	volumes = sp.zeros((numberCombos, 3))

	print "Generating volumes for: ", numberCombos, "of combinations"

	def index_coordinate_assignment(index_coordinate):
		volumes[index_coordinate[0]] = MotifsLeastSquares.findVolumeByCoordinates(index_coordinate[1])
		if (index_coordinate[0] % 1000 == 0): print "Generated ", index_coordinate[0], "volumes"

	pool = multiprocessing.Pool(processes=4)
	#pool.map(index_coordinate_assignment, enumerate(aminoAcidCoordinateMotifs(coordList)))
	volumes = sp.array(pool.map( MotifsLeastSquares.findVolumeByCoordinates, aminoAcidCoordinateMotifs(coordList)))
	#map(index_coordinate_assignment, enumerate(aminoAcidCoordinateMotifs(coordList)))
	return volumes

if __name__ == "__main__":
	coordList_native = PdbUtilityFunctions.extractCoordinateListFromPdb(PDB_FILE_NATIVE)
	coordList_generated = PdbUtilityFunctions.extractCoordinateListFromPdb(PDB_FILE_GENERATED)

	volume_native = volumeOfMotifs(coordList_native)
	volume_generated = volumeOfMotifs(coordList_generated)