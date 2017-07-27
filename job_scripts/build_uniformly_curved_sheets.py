#!/usr/bin/env python2.7
import os
import sys

import numpy as np

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

from rosetta.protocols.parametric_design import *


def make_sheets(data_path):
    pass

def make_a_sheet(sheet_path, num_strands, strand_length, crease_strand_angle, folding_angle):
    '''Make a beta sheet and record the statistics of this sheet.'''
   
    ## Build the sheet

    skeleton = make_uniformly_curved_rectangular_antiparallel_skeleton(num_strands,
                strand_length, crease_strand_angle, folding_angle)

    psd = ParametricSheetDesigner()
    psd.skeleton(skeleton)
    psd.set_constraint_weight(1)

    sheet = rosetta.core.pose.Pose()
    psd.apply(sheet)

    sheet.dump_pdb(os.path.join(sheet_path, 'sheet.pdb'))
   
    ## Evaluate the quality of the sheet 
 
    evaluate_sheet_quality(sheet_path, sheet, num_strands, strand_length, skeleton)

def evaluate_sheet_quality(sheet_path, sheet, num_strands, strand_length, skeleton):
    '''Save the quality metrics of a sheet into json files.'''

    # Get the strand sequence positions
    
    strands = [(i * strand_length + 1, (i + 1) * strand_length) for i in range(num_strands)]

    # Create the BetaSheetQualityFilter
    		
    bsqf = filters.BetaSheetQualityFilter()
    
    for s in strands:
        bsqf.add_a_strand(s[0], s[1])
    
    # Calculate the N-CA-C angle deviations

    deviations = bsqf.n_ca_c_angle_deviations(sheet);
    deviations_degree = [np.degrees(x) for x in deviations]

    # Save the N-CA-C angle deviations into a TSV file 
    
    ll = [[i + 1, deviations_degree[i]] for i in range(len(deviations_degree))]
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'deviations_degree.tsv'), sorted(ll, key=lambda x : x[1]))

    # Calculate the Ramachandran scores
    
    rama_scores = list(bsqf.rama_scores(sheet))

    # Save the Ramachandran scores into a JSON file 
    
    ll = [[i + 1, rama_scores[i]] for i in range(len(rama_scores))]
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'rama_scores.tsv'), sorted(ll, key=lambda x : x[1]))

    # Calculate distances of the nearest neighbors
    
    nearest_neighbor_distances = list(bsqf.nearest_neighbor_distances(sheet))

    # Save the distances of nearst neighbors into a JSON file 
    
    ll = [[i + 1, nearest_neighbor_distances[i]] for i in range(len(nearest_neighbor_distances))]
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'nearest_neighbor_distances.tsv'), sorted(ll, key=lambda x : x[1]))
    
    # Calculate energies of expected H-bonds
    
    sheet_hbonds = skeleton.topology().intra_sheet_h_bonds()
    bsqf.set_hbonds(sheet_hbonds)
    hbond_scores = bsqf.hbond_scores(sheet)

    # Save energies of expected H-bonds into a JSON file

    ll = []
    for i in range(1, len(sheet_hbonds) + 1):
        ll.append([int(sheet_hbonds[i].donor_res), int(sheet_hbonds[i].acceptor_res), hbond_scores[i]])
    
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'hbond_scores.tsv'), sorted(ll, key=lambda x : x[2]))

if __name__ == '__main__':

    data_path = sys.argv[1]
    input_path = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1] 
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 4:
        num_jobs = int(sys.argv[3])
        job_id = int(sys.argv[4]) - 1

    pyrosetta.init()

    make_a_sheet(data_path, 3, 7, np.radians(30.0), np.radians(50.0))
