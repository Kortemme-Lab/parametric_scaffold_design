#!/usr/bin/env python2.7
import sys
import os

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

    for mj in movable_jumps:
        ssa.add_movable_jump(mj)

    for c in connections:
        ssa.add_connection(seqpos_map[c[0]], seqpos_map[c[1]], c[2])

    ssa.apply(pose)


if __name__ == '__main__':

    pyrosetta.init()

    #save_residue_transformation('data/test_assemble/after_fast_design.pdb', ('C', 34), ('A', 11), 'data/test_assemble/T_C34_A11.json')
   
    pdb_file1 = 'data/antiparallel_sheets/4_0_50_10/sheet.pdb'
    pdb_file2 = 'data/straight_helices/10/helix.pdb'

    pose1 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose1, pdb_file1)
    pose2 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose2, pdb_file2)

    seqpos_map = make_seqpos_map(pose1, pose2)

    T = PPSD.io.load_rigid_body_transformation_from_file('data/test_assemble/T_C34_A11.json')

    pre_assemble(pose1, pose2, ('B', 10), ('A', 5), T)
    
    pose1.dump_pdb('data/test_assemble/after_pre_assemble.pdb')

    connections = [((1, 'A', 7), (1, 'B', 8), 4),
                   ((1, 'B', 14), (1,'C', 15), 4),
                   ((2, 'A', 10), (1,'A', 1), 4)]

    assemble(pose1, [3], connections, seqpos_map)

    pose1.dump_pdb('data/test_assemble/after_assemble.pdb')
