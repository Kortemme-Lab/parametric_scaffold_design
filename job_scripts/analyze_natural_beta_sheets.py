#!/usr/bin/env python2.7

import numpy as np

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

from rosetta.protocols.parametric_design import *


def extract_natural_sheet_info(pdb_path):
    '''Extract statistics from a natural beta sheet 
    stored in a PDB file
    '''
    # Read the sheet from PDB

    sheet = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(sheet, pdb_path);

    # Create the BetaSheetQualityFilter
    		
    bsqf = filters.BetaSheetQualityFilter()

    # Add the whole pose as a single strand

    bsqf.add_a_strand(1, sheet.size())

    # Calculate the N-CA-C angle deviations

    deviations = bsqf.n_ca_c_angle_deviations(sheet);
    deviations_degree = [np.degrees(x) for x in deviations]

    # Calculate the Ramachandran scores
    rama_scores = list(bsqf.rama_scores(sheet))

    # Calculate the distances of nearest neighbors
    nearest_neighbor_distances = list(bsqf.nearest_neighbor_distances(sheet))

    # Calculate scores of all H-bonds
    all_hb_scores = list(bsqf.all_hbond_scores(sheet))

    return deviations_degree, rama_scores, nearest_neighbor_distances, all_hb_scores


if __name__ == '__main__':

    pyrosetta.init()

    pdb_path = 'test_sheet.pdb'
    deviations_degree, rama_scores, nearest_neighbor_distances, all_hb_scores = \
        extract_natural_sheet_info(pdb_path)

    PPSD.plot.plot_histogram(deviations_degree, 'N-CA-C angle deviations (degree)')
