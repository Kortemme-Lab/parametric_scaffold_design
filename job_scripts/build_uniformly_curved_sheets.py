#!/usr/bin/env python2.7
import os
import sys
import json

import numpy as np
import matplotlib.pyplot as plt

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

from rosetta.protocols.parametric_design import *


def make_sheets(data_path, strand_directions, strand_length, crease_strand_angles_degree, folding_angles_degree):
    '''Make multiple beta sheets and save each into a separate directory.'''
    for i, csa in enumerate(sorted(crease_strand_angles_degree)):
        for j, fa in enumerate(sorted(folding_angles_degree)):
            sheet_path = os.path.join(data_path, '_'.join([str(i), str(j), str(int(csa)), str(int(fa))]))
            make_a_sheet(sheet_path, strand_directions, strand_length, np.radians(csa), np.radians(fa))


def make_a_sheet(sheet_path, strand_directions, strand_length, crease_strand_angle, folding_angle):
    '''Make a beta sheet and record the statistics of this sheet.'''
    ## Create directory
    if not os.path.exists(sheet_path):
        os.mkdir(sheet_path)

    ## Build the sheet

    skeleton = make_uniformly_curved_rectangular_skeleton(strand_directions,
                strand_length, crease_strand_angle, folding_angle)

    psd = ParametricSheetDesigner()
    psd.skeleton(skeleton)
    psd.set_constraint_weight(1)

    sheet = rosetta.core.pose.Pose()
    psd.apply(sheet)

    sheet.dump_pdb(os.path.join(sheet_path, 'sheet.pdb'))
   
    ## Evaluate the quality of the sheet 
 
    evaluate_sheet_quality(sheet_path, sheet, len(strand_directions), strand_length, skeleton)

    # Save the information of the sheet into a JSON file

    with open(os.path.join(sheet_path, 'info.json'), 'w') as f:
        json.dump({'strand_directions': list(strand_directions),
                   'strand_length' : strand_length,
                   'crease_strand_angles_degree' : np.degrees(crease_strand_angle),
                   'folding_angles_degree' : np.degrees(folding_angle)}, f)
    
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

    # Save the Ramachandran scores into a TSV file 
    
    ll = [[i + 1, rama_scores[i]] for i in range(len(rama_scores))]
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'rama_scores.tsv'), sorted(ll, key=lambda x : x[1]))

    # Calculate distances of the nearest neighbors
    
    nearest_neighbor_distances = list(bsqf.nearest_neighbor_distances(sheet))

    # Save the distances of nearest neighbors into a TSV file 
    
    ll = [[i + 1, nearest_neighbor_distances[i]] for i in range(len(nearest_neighbor_distances))]
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'nearest_neighbor_distances.tsv'), sorted(ll, key=lambda x : x[1]))
    
    # Calculate energies of expected H-bonds
    
    sheet_hbonds = skeleton.topology().intra_sheet_h_bonds()
    bsqf.set_hbonds(sheet_hbonds)
    hbond_scores = bsqf.hbond_scores(sheet)

    # Save energies of expected H-bonds into a TSV file

    ll = []
    for i in range(1, len(sheet_hbonds) + 1):
        ll.append([int(sheet_hbonds[i].donor_res), int(sheet_hbonds[i].acceptor_res), hbond_scores[i]])
    
    PPSD.io.ll_to_tsv(os.path.join(sheet_path, 'hbond_scores.tsv'), sorted(ll, key=lambda x : x[2]))

def read_sheet_info_into_a_dict(sheet_path):
    '''Read all information and evaluation about a sheet into a dictionary.'''
    info_dict = {}
    
    # Load the basic information
    
    with open(os.path.join(sheet_path, 'info.json'), 'r') as f:
        info_dict = json.load(f)

    # Load the N-CA-C angle deviations 
    
    info_dict['deviations_degree'] = PPSD.io.tsv_to_ll(os.path.join(sheet_path, 'deviations_degree.tsv'), [lambda x : int(x), lambda x : float(x)])

    # Load Ramachandran scores
    
    info_dict['rama_scores'] = PPSD.io.tsv_to_ll(os.path.join(sheet_path, 'rama_scores.tsv'), [lambda x : int(x), lambda x : float(x)])

    # Load the nearest neighbor distances 
    
    info_dict['nearest_neighbor_distances'] = PPSD.io.tsv_to_ll(os.path.join(sheet_path, 'nearest_neighbor_distances.tsv'), [lambda x : int(x), lambda x : float(x)])

    # Load the H-bond scores
    
    info_dict['hbond_scores'] = PPSD.io.tsv_to_ll(os.path.join(sheet_path, 'hbond_scores.tsv'), [lambda x : int(x), lambda x : int(x), lambda x : float(x)])

    return info_dict

def make_summary_figures(data_path):
    '''Make figures that summarize the set of built sheets.''' 

    def x(dn):
        '''Get the x index of a sheet given its directory name'''
        return int(dn.split('_')[0])

    def y(dn):
        '''Get the y index of a sheet given its directory name'''
        return int(dn.split('_')[1])

    num_x = max([x(dn) for dn in os.listdir(data_path) if dn not in ['figures']]) + 1
    num_y = max([y(dn) for dn in os.listdir(data_path) if dn not in ['figures']]) + 1

    # Save the information dictionaries into a 2D array

    info_dicts = [[None for i in range(num_y)] for j in range(num_x)]

    for s in os.listdir(data_path):
        if s not in ['figures']:
            info_dicts[x(s)][y(s)] = read_sheet_info_into_a_dict(os.path.join(data_path, s))

    # Make a directory for figures

    figure_path = os.path.join(data_path, 'figures')
    if not os.path.exists(figure_path):
        os.mkdir(figure_path)

    # Plot the largesest N-CA-C angle deviations

    make_heat_map(info_dicts, lambda d : max(x[1] for x in d['deviations_degree']), 
            'largest N-CA-C angle deviations', figure_path)

    # Plot the smallest N-CA-C angle deviations

    make_heat_map(info_dicts, lambda d : min(x[1] for x in d['deviations_degree']), 
            'smallest N-CA-C angle deviations', figure_path)

    # Plot the largesest nearest neighbor distances

    make_heat_map(info_dicts, lambda d : max(x[1] for x in d['nearest_neighbor_distances']), 
            'largest nearest neighbor distances', figure_path)

    # Plot the smallest nearest neighbor distances

    make_heat_map(info_dicts, lambda d : min(x[1] for x in d['nearest_neighbor_distances']), 
            'smallest nearest neighbor distances', figure_path)

    # Plot the largesest Ramachandran scores

    make_heat_map(info_dicts, lambda d : max(x[1] for x in d['rama_scores']), 
            'largest rama scores', figure_path)

    # Plot the largest H-bond scores

    make_heat_map(info_dicts, lambda d : max(x[2] for x in d['hbond_scores']), 
            'largest hbond scores', figure_path)

def make_heat_map(info_dicts, get_data, title='', fig_save_path=None):
    '''Make a heat map to show certain metric of the built sheets.
    The get_data argument is a function to get specific data from an
    information dictionary.
    '''
    num_x = len(info_dicts)
    num_y = len(info_dicts[0])

    data = [[get_data(d) for d in l ] for l in info_dicts]
 
    plt.close()

    fig, ax = plt.subplots()

    cax = ax.imshow(np.transpose(data), origin='lower', interpolation='nearest')
    cbar = fig.colorbar(cax)

    plt.xticks(range(num_x), [info_dicts[i][0]['crease_strand_angles_degree'] for i in range(num_x)])
    plt.yticks(range(num_y), [info_dicts[0][i]['folding_angles_degree'] for i in range(num_y)])
    
    plt.xlabel('Crease-strand angle (degree)')
    plt.ylabel('folding angle (degree)')

    plt.title(title)

    if fig_save_path:
        plt.savefig(os.path.join(fig_save_path, '_'.join(title.split()) + '.png' ))
    else:
        plt.show()

    

if __name__ == '__main__':

    data_path = sys.argv[1]
    input_path = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1] 
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 4:
        num_jobs = int(sys.argv[3])
        job_id = int(sys.argv[4]) - 1

    pyrosetta.init()

    strand_directions = rosetta.utility.vector1_bool()
    strand_directions.append(True)
    strand_directions.append(True)
    strand_directions.append(True)

    make_sheets(data_path, strand_directions, 7, [10, 20, 30, 40, 50], [10, 20, 30, 40, 50])

    make_summary_figures(data_path)
