#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Constants
import MotifsWithLowestStd
import MotifsLeastSquares
import PdbUtilityFunctions

import scipy as sp
import random
import cPickle

PDB_FILE = "data/3ZOB-one.pdb"
OUTPUT_FILE = "hamidMethod/hamidMethod-3ZOB-gen-tub.pickle"

def printStats(stats):
    i = 1
    for coord in stats:
        print "Residue %d:" % i
        print u"\tx: μ={0:.2f}, σ={1:.2f}".format(coord[0]['mean'], coord[0]['std'])
        print u"\ty: μ={0:.2f}, σ={1:.2f}".format(coord[1]['mean'], coord[1]['std'])
        print u"\tz: μ={0:.2f}, σ={1:.2f}".format(coord[2]['mean'], coord[2]['std'])
        i += 1

def hamidMethod():
    # Preliminary step: We have to look at the PDB file and identify the
    # motifs that are most hydrophobic from the protein. (Find M)
    n_size_of_protein = PdbUtilityFunctions.findLengthOfProteinFromPdb(PDB_FILE)
    m_hydrophobic = PdbUtilityFunctions.findHydrophobicAcidsNum(PDB_FILE)

    print "Analyzing:", PDB_FILE, "with n = %d amino acids, of which m = %d are hydrophobic" \
                                  % (n_size_of_protein, m_hydrophobic)

    # First, create an initial model using our least squares model, as we've
    # done before.
    list_of_motifs_with_volumes = MotifsWithLowestStd.findListOfMotifsWithVolumes(PDB_FILE)
    filtered_output = MotifsWithLowestStd.filterOutUnneededMotifs(
        list_of_motifs_with_volumes, n_size_of_protein)
    print "Visited AA indices: ", filtered_output["visited"]
    initial_coordinates = MotifsLeastSquares.findCoordinates(filtered_output["motifs"],
                                                             filtered_output["visited"])

    # Store the initial coordinates.
    current_coordinates = initial_coordinates
    prev_motifs = filtered_output["motifs"]
    coordinates = []

    # In a loop up to 20 times:
    for i in xrange(100):
        print "Current Iteration: ", i+1

        # Next, we have to a quick analysis to check which of the least squares
        # motifs fit the most with the model (3M - 6).  Choose 3N - 3M motifs, randomly
        # selected from the *protein* model, take their volumes, and randomly change
        # the values +-10%.
        print "Sorting the motifs in the model by difference with volume mean of population..."
        motif_volume_difference = []
        new_motifs = []
        for motif in prev_motifs:
            # Find the volume of the volume
            try:
                motifs_coordinates = MotifsLeastSquares.getCoordinatesByAminoAcidIndex(
                    motif["index"], current_coordinates, filtered_output["visited"]
                )
            except Exception:
                import pdb; pdb.set_trace()

            vol_model = MotifsLeastSquares.findVolumeByCoordinates(motifs_coordinates)
            motif_volume_difference.append({
                "index": motif["index"], "modelV": vol_model, "diff": abs(vol_model - motif["unVol"])
            })
        motif_volume_difference = sorted(motif_volume_difference, key=lambda x: x["diff"])

        print "Appending the bottom 3m - 6 = %d to our new motif set..." % (3*m_hydrophobic - 6)
        for i in xrange(3*m_hydrophobic - 6):
            motif = motif_volume_difference[i]
            new_motifs.append({"unVol": motif["modelV"], "index": motif["index"]})

        print "Choosing 3n - 3m = %d random motifs from the model..." % (3*n_size_of_protein - 3*m_hydrophobic)
        random_motif_indices = []
        for i in xrange(3*n_size_of_protein - 3*m_hydrophobic):
            random_motif_indices.append(random.sample(filtered_output["visited"], 4))

        print "Changing the value of each by +- 10% and appending those motifs to the motif set..."
        for motif_index in random_motif_indices:
            motifs_coordinates = MotifsLeastSquares.getCoordinatesByAminoAcidIndex(
                motif_index, current_coordinates, filtered_output["visited"]
            )
            vol_model = MotifsLeastSquares.findVolumeByCoordinates(motifs_coordinates)
            # Add +- 10% to the volume
            vol_model_r = vol_model * random.uniform(.9, 1.1)
            new_motifs.append({
                "index": motif["index"], "unVol": vol_model_r
            })

        # Perform least squares.
        current_coordinates = MotifsLeastSquares.findCoordinatesFromInitial(new_motifs,
                                                                 filtered_output["visited"],
                                                                 current_coordinates)

        # Chunk it into coordinates
        current_coordinates_grouped = []
        for i in range(0, len(current_coordinates), 3):
            current_coordinates_grouped.append((current_coordinates[i], current_coordinates[i + 1],
                                        current_coordinates[i + 2]))

        # Store the coordinates
        coordinates.append(current_coordinates_grouped)
        prev_motifs = new_motifs

        print ""

    # Do processing on the std. deviation and average of the coordinates
    print "Processing stats on coordinates..."
    stat_container = []
    for i in xrange(n_size_of_protein):
        stat_container.append(([], [], []))

    for coord_iter in coordinates:
        for coord_index in xrange(len(coord_iter)):
            stat_container[coord_index][0].append(coord_iter[coord_index][0])
            stat_container[coord_index][1].append(coord_iter[coord_index][1])
            stat_container[coord_index][2].append(coord_iter[coord_index][2])

    stats = []
    for coord in stat_container:
        stats.append(({"mean": sp.mean(coord[0]), "std": sp.std(coord[0])},
                      {"mean": sp.mean(coord[1]), "std": sp.std(coord[1])},
                      {"mean": sp.mean(coord[2]), "std": sp.std(coord[2])}))

    printStats(stats)

    return {"stats": stats, "initial": initial_coordinates, "iter_coords": coordinates,
            "visited": filtered_output["visited"] }

if __name__ == "__main__":
    returnObj = hamidMethod()

    with open(OUTPUT_FILE, 'w') as f:
        cPickle.dump(returnObj, f)

