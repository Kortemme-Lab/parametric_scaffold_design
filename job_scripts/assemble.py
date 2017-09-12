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


def pre_assemble(pdb_file1, pdb_file2, res1, res2, frame_transformation):
    '''Assemble two PDB files into one file given the relative
    orientation of one residue from each structures.
    
    Residues are specified by a tuple of (chain_id, residue_id).
    
    Return the pre-assembled pose.
    '''
    
    pose1 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose1, pdb_file1)
    pose2 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose2, pdb_file2)

    seqpos1 = pose1.pdb_info().pdb2pose(res1[0], res1[1])
    seqpos2 = pose2.pdb_info().pdb2pose(res2[0], res2[1]) + pose1.size()
    
    #print_pose_pdb_map(pose1)###DEBUG 
    #print_pose_pdb_map(pose2)###DEBUG

    # Merge pose2 into pose1

    pose1.append_pose_by_jump(pose2, pose1.size())
    
    #print_pose_pdb_map(pose1)###DEBUG

    transform_chain(pose1, seqpos1, seqpos2, frame_transformation)

    pose1.dump_pdb('data/test_assemble/test.pdb')


if __name__ == '__main__':

    pyrosetta.init()

    #save_residue_transformation('data/test_assemble/after_fast_design.pdb', ('C', 34), ('A', 11), 'data/test_assemble/T_C34_A11.json')

    T = PPSD.io.load_rigid_body_transformation_from_file('data/test_assemble/T_C34_A11.json')

    pre_assemble('data/antiparallel_sheets/4_4_50_50/sheet.pdb', 'data/straight_helices/10/helix.pdb', ('B', 10), ('A', 5), T)
