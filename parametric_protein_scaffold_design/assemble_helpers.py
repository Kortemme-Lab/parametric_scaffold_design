import pyrosetta
from pyrosetta import rosetta
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

def loops_to_rosetta_loops(loops):
    '''Create a rosetta Loops object from a list of pairs (start, stop).'''
    loops_rosetta = rosetta.protocols.loops.Loops()
    for l in loops:
        loops_rosetta.add_loop(rosetta.protocols.loops.Loop(l[0], l[1], l[1], 0, True))

    return loops_rosetta

def fast_loop_build(pose, loops):
    '''Build loops quickly. Loops are given as pairs of (start, stop)'''
    loops_rosetta = loops_to_rosetta_loops(loops)

    loop_modeler = rosetta.protocols.loop_modeler.LoopModeler()
    loop_modeler.set_loops(loops_rosetta)
    loop_modeler.setup_kic_config()
    loop_modeler.fullatom_stage().mark_as_test_run()
    loop_modeler.centroid_stage().mark_as_test_run()

    loop_modeler.apply(pose)
    return loop_modeler.build_stage().was_successful()

def minimization_loop_closing(pose, loops):
    '''Close the loops of a pose by minimization.'''
    pose.dump_pdb('debug/before_minimize.pdb')
    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(1, 14, 1)
    ft.add_edge(1, 38, 2)
    ft.add_edge(1, 48, 3)
    ft.add_edge(1, 12, -1)
    ft.add_edge(14, 13, -1)
    ft.add_edge(14, 36, -1)
    ft.add_edge(38, 37, -1)
    ft.add_edge(38, 46, -1)
    ft.add_edge(48, 47, -1)
    ft.add_edge(48, 54, -1)
    pose.fold_tree(ft)
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    min_mover = rosetta.protocols.simple_moves.MinMover()
    mm = rosetta.core.kinematics.MoveMap()
    mm.set_jump(1, True)

    #mm.set(rosetta.core.id.THETA, True)
    #mm.set(rosetta.core.id.D, True)
    for loop in loops:
        for seqpos in range(loop[0], loop[1] + 1):
            #print pose.fold_tree().is_cutpoint( seqpos )
            mm.set_bb(seqpos, True)
    
    min_opts = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True )

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn.set_weight(rosetta.core.scoring.linear_chainbreak, 10)
    #sfxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015_cart")


    #min_mover.cartesian(True);
    min_mover.score_function(sfxn)
    min_mover.movemap(mm)
    min_mover.min_options(min_opts)

    min_mover.apply(pose)

    fast_loop_build(pose, loops)

def mutate_residues(pose, res_list, aa_list):
    '''Mutate a list of residues. The list of AAs could
    either be 1 letter code or 3 letter code.
    '''
    aa_name_map = {'A':'ALA', 'P':'PRO', 'V':'VAL', 'L':'LEU', 'I':'ILE', 'M':'MET',
                   'F':'PHE', 'Y':'TYR', 'W':'TRP', 'S':'SER', 'T':'THR', 'C':'CYS',
                   'K':'LYS', 'R':'ARG', 'H':'HIS', 'D':'ASP', 'E':'GLU', 'N':'ASN',
                   'Q':'GLN', 'G':'GLY'}


    mutater = rosetta.protocols.simple_moves.MutateResidue()
    for i in range(len(res_list)):
        name = aa_list[i] if len(aa_list[i]) == 3 else aa_name_map[aa_list[i]]
        mutater.set_res_name(name)
        mutater.set_target(res_list[i])
        mutater.apply(pose)
