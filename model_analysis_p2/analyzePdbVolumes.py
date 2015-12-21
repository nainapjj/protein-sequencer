import itertools
import math
import multiprocessing

import scipy as sp
import scipy.optimize
import scipy.stats

from ..model_generation_p2 import Constants, MotifsLeastSquares, MotifsWithLowestStd, StuartFindAllPossibleMotifs
from ..model_util_p2 import PdbUtilityFunctions

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
    volumes = sp.array(pool.map(MotifsLeastSquares.findVolumeByCoordinates, aminoAcidCoordinateMotifs(coordList)))
    #map(index_coordinate_assignment, enumerate(aminoAcidCoordinateMotifs(coordList)))
    return volumes


def filterInOnlyAllHydrophobic(motif):
    hydrophobic = PdbUtilityFunctions.returnHydrophobicAminos()
    hydrophobicMap = StuartFindAllPossibleMotifs.dict
    singleLetterHydro = [hydrophobicMap[hydro.strip()] for hydro in hydrophobic]

    for aa in motif["amino"]:
        if aa not in singleLetterHydro:
            return False

    return True


def volumeErrorsForMotifs(pdbModelFile, filterFunc=None, output=False, list_of_motifs_with_volumes=None):
    # Extract the residual sequence of the protein
    resSequence = PdbUtilityFunctions.getResidualSequence(pdbModelFile)

    # Get the size of the protein
    n_size_of_protein = PdbUtilityFunctions.findLengthOfProteinFromPdb(pdbModelFile)

    # Extract the alpha carbon coordinates from the pdb
    pdbCoordinates = PdbUtilityFunctions.extractCoordinateListFromPdb(pdbModelFile)

    # Next, obtain the motifs we want to analyze
    if list_of_motifs_with_volumes is None:
        list_of_motifs_with_volumes = MotifsWithLowestStd.findListOfMotifsWithVolumes(pdbModelFile)
    filtered_motifs = MotifsWithLowestStd.filterOutUnneededMotifs(
        list_of_motifs_with_volumes, n_size_of_protein)

    # Now perform the filter step
    if filterFunc:
        filtered_motifs = filter(filterFunc, filtered_motifs["motifs"])
    else:
        filtered_motifs = filtered_motifs["motifs"]

    # Now we could perform the volume difference analysis
    deltas = []
    stds = []
    actualVolumes = []
    expectedVolumes = []
    indices = []
    coordinates = []

    for motif in filtered_motifs:
        # Convert the pdb-based index to a zero-based index
        mot1Ind = PdbUtilityFunctions.fromResIndToZeroInd(motif["index"][0], resSequence)
        mot2Ind = PdbUtilityFunctions.fromResIndToZeroInd(motif["index"][1], resSequence)
        mot3Ind = PdbUtilityFunctions.fromResIndToZeroInd(motif["index"][2], resSequence)
        mot4Ind = PdbUtilityFunctions.fromResIndToZeroInd(motif["index"][3], resSequence)

        try:
            # Find the volume of the motif
            volArea = MotifsLeastSquares.findVolumeByCoordinates((pdbCoordinates[mot1Ind], pdbCoordinates[mot2Ind],
                                                                  pdbCoordinates[mot3Ind], pdbCoordinates[mot4Ind]))
        except IndexError:
            import pdb; pdb.set_trace()

        # Compare it to the mean volume, after indexating
        delta = volArea - motif["unVol"]

        # Append the data files
        actualVolumes.append(volArea)
        expectedVolumes.append(motif["unVol"])
        deltas.append(delta)
        stds.append(motif["std"])
        indices.append(motif["index"])
        coordinates.append((pdbCoordinates[mot1Ind], pdbCoordinates[mot2Ind],
                            pdbCoordinates[mot3Ind], pdbCoordinates[mot4Ind]))

        if output:
            # Output the result
            print "motif: %s, std: %.2f actual volume: %.2f, expected volume: %.2f, " \
                  "error: %.2f" % (motif["amino"], motif["std"], volArea, motif["unVol"], delta)

    return {"aVol": actualVolumes, "eVol": expectedVolumes, "deltas": deltas, "stds": stds,
            "indices": indices, "coordinates": coordinates}


def deleteByIndex(dict, indexes):
    for key in dict.keys():
        for index in sorted(indexes, reverse=True):
            del(dict[key][index])


def maximizeCorrelation():

    motifs_for_3ZOB = MotifsWithLowestStd.findListOfMotifsWithVolumes("data/3ZOB-one.pdb")

    def correlationFitness((mean_pairwise, min_pairwise, min_actual_vol, max_actual_vol)):
        print "Evaluating for: ", mean_pairwise, min_pairwise, min_actual_vol, max_actual_vol
        Constants.MEAN_PAIRWISE = mean_pairwise
        Constants.MIN_PAIRWISE = min_pairwise

        data_3ZOB = volumeErrorsForMotifs("data/3ZOB-one.pdb", list_of_motifs_with_volumes=motifs_for_3ZOB)
        #data_1DF4 = volumeErrorsForMotifs("pdbs/1DF4.pdb")
        #data_2RJV = volumeErrorsForMotifs("pdbs/2RJV.pdb")
        #data_2RJY = volumeErrorsForMotifs("pdbs/2RJY.pdb")

        # Filter out the motifs with unusual actual volumes
        deleteByIndex(data_3ZOB, sp.where((sp.array(data_3ZOB["aVol"]) <= min_actual_vol) |
                                          (sp.array(data_3ZOB["aVol"]) >= max_actual_vol))[0])

        pearson_3ZOB = scipy.stats.pearsonr(data_3ZOB["aVol"], data_3ZOB["eVol"])[0]

        if len(data_3ZOB["aVol"]) < 5 or math.isnan(pearson_3ZOB):
            pearson_3ZOB = 0

        print "Pearson: %f" % pearson_3ZOB
        return -pearson_3ZOB

    bounds = [(0, 20), (0, 20), (4, 10), (10, 30)]
    minimizer_kwargs = {"method": "L-BFGS-B", "bounds": bounds}
    result = scipy.optimize.basinhopping(correlationFitness, [Constants.MEAN_PAIRWISE, Constants.MIN_PAIRWISE,
                                                              Constants.MIN_ACTUAL_VOL, Constants.MAX_ACTUAL_VOL],
                                         niter=100, disp=True, minimizer_kwargs=minimizer_kwargs, stepsize=3)
    print result

    return result



def main_test():
    print "Original PDB"
    return volumeErrorsForMotifs('pdbs/1DF4.pdb')

def main_compare():
    #volumeErrorsForMotifs("tinker/1DF4-gen-bb-tinker-mini.pdb")
    print "Data for 2RJY"
    print ""

    print "Original PDB"
    volumeErrorsForMotifs('pdbs/2RJY.pdb')

    print ""
    print "Generated PDB"
    volumeErrorsForMotifs('hamidMethod/2RJY-gen-tub.pdb')

    print ""
    print "Tinker PDB after Minimization"
    volumeErrorsForMotifs('tinker/2RJY-gen-tub-tinker-mini.pdb')

def main_analyzeVolume():
    coordList_native = PdbUtilityFunctions.extractCoordinateListFromPdb(PDB_FILE_NATIVE)
    coordList_generated = PdbUtilityFunctions.extractCoordinateListFromPdb(PDB_FILE_GENERATED)

    volume_native = volumeOfMotifs(coordList_native)
    volume_generated = volumeOfMotifs(coordList_generated)

if __name__ == "__main__":
    maximizeCorrelation()