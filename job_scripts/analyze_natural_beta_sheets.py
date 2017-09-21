#!/usr/bin/env python2.7

import os
import sys

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

def save_one_natural_sheet_info(data_path, pdb_path):
    '''Save the statistics of one beta sheet into data files.'''
    # Extract the info

    try:
        deviations_degree, rama_scores, nearest_neighbor_distances, all_hb_scores = \
            extract_natural_sheet_info(pdb_path)
    except:
        return

    # Save the info into different files 

    if len(deviations_degree) > 0:
        PPSD.io.safe_append(os.path.join(data_path, 'deviations_degree.txt'),
                '\n'.join(["{0:0.2f}".format(x) for x in deviations_degree]) + '\n')

    if len(rama_scores) > 0:
        PPSD.io.safe_append(os.path.join(data_path, 'rama_scores.txt'),
                '\n'.join(["{0:0.2f}".format(x) for x in rama_scores]) + '\n')

    if len(nearest_neighbor_distances) > 0:
        PPSD.io.safe_append(os.path.join(data_path, 'nearest_neighbor_distances.txt'),
                '\n'.join(["{0:0.2f}".format(x) for x in nearest_neighbor_distances]) + '\n')

    if len(all_hb_scores) > 0:
        PPSD.io.safe_append(os.path.join(data_path, 'all_hb_scores.txt'),
                '\n'.join(["{0:0.2f}".format(x) for x in all_hb_scores]) + '\n')

def save_info_of_natrual_sheets_from_cath(data_path, sheets_path, num_jobs, job_id):
    '''Save information of natural sheets from the cath database.'''
    my_super_families = PPSD.job_distributors.SGEJobDistributor.get_sub_list(
            os.listdir(sheets_path), num_jobs, job_id)

    for sf in my_super_families:
        for p in os.listdir(os.path.join(sheets_path, sf)):
            save_one_natural_sheet_info(data_path, os.path.join(sheets_path, sf, p))

def plot_natural_sheet_info(data_path, savefig=False):
    '''Plot the statistics of natural beta sheets.'''
    
    with open(os.path.join(data_path, 'deviations_degree.txt'), 'r') as f:
        PPSD.plot.plot_histogram([float(x) for x in f.read().split()], 'N-CA-C angle deviations (degree)',
                x_min=-15, x_max=15, savefig_path=data_path if savefig else None)

    with open(os.path.join(data_path, 'rama_scores.txt'), 'r') as f:
        PPSD.plot.plot_histogram([float(x) for x in f.read().split()], 'Ramachandran scores',
                x_min=0, x_max=15, savefig_path=data_path if savefig else None)

    with open(os.path.join(data_path, 'nearest_neighbor_distances.txt'), 'r') as f:
        PPSD.plot.plot_histogram([float(x) for x in f.read().split()], 'Distance to nearst neighbor',
                x_min=0, x_max=10, savefig_path=data_path if savefig else None)

    with open(os.path.join(data_path, 'all_hb_scores.txt'), 'r') as f:
        PPSD.plot.plot_histogram([float(x) for x in f.read().split()], 'H-Bond scores', 
                x_min=-2, x_max=0, savefig_path=data_path if savefig else None)


if __name__ == '__main__':

    data_path = sys.argv[1]
    input_path = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1] 
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 4:
        num_jobs = int(sys.argv[3])
        job_id = int(sys.argv[4]) - 1

    #pyrosetta.init()
    #save_info_of_natrual_sheets_from_cath(data_path, input_path, num_jobs, job_id)

    plot_natural_sheet_info(data_path, True)
