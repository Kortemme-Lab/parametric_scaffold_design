from rosetta.protocols.assembly.secondary_structure_assembly import *


def save_residue_transformation(pdb_file, res1, res2, fout):
    '''Save the rigid body transformation from the res1 frame to the res2 frame.
    Residues are specified by a tuple of (chain_id, residue_id).
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    
    seqpos1 = pose.pdb_info().pdb2pose(res1[0], res1[1])
    seqpos2 = pose.pdb_info().pdb2pose(res2[0], res2[1])

    T = calc_residue_frame_transformation(pose, seqpos1, seqpos2)
    IO.write_rigid_body_transformation_to_file(T, fout)

def make_seqpos_map(poses):
    '''Make a map from the (structure_id, chain, pdb_number)
    to the merged pose. The structure_id starts from 0.
    '''
    seqpos_map = {}

    infos = [pose.pdb_info() for pose in poses]

    offset = 0
    for i in range(len(poses)):
        for j in range(1, poses[i].size() + 1):
            seqpos_map[(i, str(infos[i].chain(j)), int(infos[i].number(j)))] = j + offset
        offset += poses[i].size()

    return seqpos_map

def pre_assemble(poses, res_pairs, frame_transformations):
    '''Assemble a list of poses into one pose given the relative
    orientation of one residue from each structures. The resulting structure
    will be saved in the first pose.
    
    Residues are specified by a pair of tuples of (pose_id, chain_id, residue_id).
    '''

    def res_to_seqpos(res, poses_before_assemble):
        '''Convert a residue into the seqpos of the assembled structure'''
        offset = 0
        for i in range(res[0]):
            offset += poses_before_assemble[i].size()

        return offset + poses_before_assemble[res[0]].pdb_info().pdb2pose(res[1], res[2])
   
    # Get the seqposes after assembly
    
    seqpos_pairs = []
    
    for res1, res2 in res_pairs:
        seqpos_pairs.append((res_to_seqpos(res1, poses), res_to_seqpos(res2, poses)))

    # Merge poses into pose1

    for i in range(1, len(poses)):
        poses[0].append_pose_by_jump(poses[i], poses[0].size())
   
    # Transform the chains

    for i in range(len(seqpos_pairs)):
        transform_chain(poses[0], seqpos_pairs[i][0], seqpos_pairs[i][1], frame_transformations[i])

