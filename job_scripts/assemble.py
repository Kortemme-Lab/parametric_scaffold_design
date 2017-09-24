#!/usr/bin/env python2.7
import sys
import os
import json
import datetime

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

from rosetta.protocols.assembly.secondary_structure_assembly import *


def print_pose_pdb_map(pose):
    '''Print the map between pose indexing and the pdb indexing.'''
    for i in range(1, pose.size() + 1):
        print i, pose.pdb_info().chain(i), pose.pdb_info().number(i)
    
    print '\n'

def save_residue_transformation(pdb_file, res1, res2, fout):
    '''Save the rigid body transformation from the res1 frame to the res2 frame.
    Residues are specified by a tuple of (chain_id, residue_id).
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    
    seqpos1 = pose.pdb_info().pdb2pose(res1[0], res1[1])
    seqpos2 = pose.pdb_info().pdb2pose(res2[0], res2[1])

    T = calc_residue_frame_transformation(pose, seqpos1, seqpos2)
    PPSD.io.write_rigid_body_transformation_to_file(T, fout)

def make_seqpos_map(pose1, pose2):
    '''Make a map from the (structure_id, chain, pdb_number)
    to the merged pose.
    '''
    seqpos_map = {}

    info1 = pose1.pdb_info()
    info2 = pose2.pdb_info()

    for i in range(1, pose1.size() + 1):
        seqpos_map[(1, str(info1.chain(i)), int(info1.number(i)))] = i

    for i in range(1, pose2.size() + 1):
        seqpos_map[(2, str(info2.chain(i)), int(info2.number(i)))] = i + pose1.size()

    return seqpos_map

def pre_assemble(pose1, pose2, res1, res2, frame_transformation):
    '''Assemble two poses into one file given the relative
    orientation of one residue from each structures.
    
    Residues are specified by a tuple of (chain_id, residue_id).
    '''
    
    seqpos1 = pose1.pdb_info().pdb2pose(res1[0], res1[1])
    seqpos2 = pose2.pdb_info().pdb2pose(res2[0], res2[1]) + pose1.size()
    
    # Merge pose2 into pose1

    pose1.append_pose_by_jump(pose2, pose1.size())
    
    transform_chain(pose1, seqpos1, seqpos2, frame_transformation)

def assemble(pose, movable_jumps, connections, seqpos_map):
    '''Assemble the secondary structures into a connected protein. 
    movable_jumps is a list of jump ids,
    connections is a list of tuples (res1, res2, connection_length) where residues are 
    specified by PDB indexing as (pose_id, chain, res_id). The PDB indices 
    are converted to Rosetta indices by the seqpos_map.
    '''
    ssa = SecondaryStructureAssembler()
    ssa.fast_design_rounds(5)

    for mj in movable_jumps:
        ssa.add_movable_jump(mj)

    for c in connections:
        ssa.add_connection(seqpos_map[c[0]], seqpos_map[c[1]], c[2])

    ssa.apply(pose)

def save_task_info(output_path, sequence, score, run_time):
    '''Save the task information into a json file.'''
    with open(os.path.join(output_path, 'task_info.json'), 'w') as f:
        json.dump({'sequence':sequence,
                   'score':score,
                   'run_time':run_time}, f)

def assemble_from_files(pdb_file1, pdb_file2, transformation_file, res1, res2, movable_jumps, connections, output_path):
    '''Assemble two secondary structures given the PDB files, transformation
    file, transformation residues, movable jumps and connections.
    The outputs will be saved to the output_path.
    '''
    start_time = datetime.datetime.now()
   
    # Load the structures

    pose1 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose1, pdb_file1)
    pose2 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose2, pdb_file2)

    seqpos_map = make_seqpos_map(pose1, pose2)

    # Do the pre-assembly

    T = PPSD.io.load_rigid_body_transformation_from_file(transformation_file)
    pre_assemble(pose1, pose2, res1, res2, T)
   
    # Do the assembly

    assemble(pose1, movable_jumps, connections, seqpos_map)

    pose1.dump_pdb(os.path.join(output_path, 'assembled.pdb'))
    PPSD.io.sequence_to_fasta_file(os.path.join(output_path, 'assembled.fasta'), 'assembled', pose1.sequence())

    end_time = datetime.datetime.now()
    run_time = end_time - start_time
    
    # Save the task info

    save_task_info(output_path, pose1.sequence(), pose1.energies().total_energy(), run_time.total_seconds())

def run_tasks(task_list, num_jobs, job_id):
    '''Run tasks that belongs to the current job thread'''
    for i, t in enumerate(task_list):
        if job_id == i % num_jobs:
            if not os.path.exists(t['output_path']): 
                os.makedirs(t['output_path'])

            assemble_from_files(t['pdb_file1'], t['pdb_file2'], t['transformation_file'], t['res1'], t['res2'],
                    t['movable_jumps'], t['connections'], t['output_path'])

def pilot_run(data_path, num_jobs, job_id):
    pdb_file1 = 'data/antiparallel_sheets/2_2_30_30/sheet.pdb'
    pdb_file2 = 'data/straight_helices/15/helix.pdb'
    transformation_file = 'database/transformations/sheet_helix_transformation.json'
    res1 = ('B', 10)
    res2 = ('A', 7)
    movable_jumps = [3]
    connections = [((1, 'A', 7), (1, 'B', 8), 2),
                   ((1, 'B', 14), (1,'C', 15), 4),
                   ((2, 'A', 15), (1,'A', 1), 4)]
    output_path = data_path

    task_list = []

    for i in range(1000):
        task_list.append({'pdb_file1':pdb_file1,
                      'pdb_file2':pdb_file2,
                      'transformation_file':transformation_file,
                      'res1':res1,
                      'res2':res2,
                      'movable_jumps':movable_jumps,
                      'connections':connections,
                      'output_path':os.path.join(output_path, str(i))})
    
    run_tasks(task_list, num_jobs, job_id)

if __name__ == '__main__':

    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    pyrosetta.init()

    #save_residue_transformation('data/test_assemble/after_fast_design.pdb', ('C', 34), ('A', 11), 'data/test_assemble/T_C34_A11.json')
   
    pilot_run(data_path, num_jobs, job_id)
