import json

import numpy as np

import pyrosetta
import pyrosetta.rosetta as rosetta

def graft_loop(pose, loop, linker_database):
    '''Graft a loop from a linker database to a pose
    at a given position. The loop is defined as (start, stop).
    '''
    # Record the target positions 
    target_ca = rosetta.numeric.xyzVector_double_t(pose.residue(loop[1]).xyz('CA'))
    target_c = rosetta.numeric.xyzVector_double_t(pose.residue(loop[1]).xyz('C'))

    # Set the loop fold tree
    
    old_ft = pose.fold_tree()
    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(1, loop[1], -1)
    ft.add_edge(1, pose.size(), 1)
    ft.add_edge(pose.size(), loop[1] + 1, -1)
    pose.fold_tree(ft)

    # Idealize the loop
    
    for i in range(loop[0], loop[1] + 1):
        rosetta.core.conformation.idealize_position(i,
                pose.conformation())
    
    with open(linker_database, 'r') as f:
        linkers = json.load(f)
    
    loop_size = loop[1] - loop[0] + 1

    for index, linker in enumerate(linkers):
        for j in range(loop_size):
                pose.set_phi(loop[0] + j, linker['phis'][j])
                if j < loop_size - 1:
                    pose.set_psi(loop[0] + j, linker['psis'][j])
                    pose.set_omega(loop[0] + j, linker['omegas'][j])

        # Check if the loop is approximately closed by the linker
    
        current_ca = pose.residue(loop[1]).xyz('CA')
        current_c = pose.residue(loop[1]).xyz('C')

        rmsd = np.sqrt(((current_ca - target_ca).length_squared() + (current_c - target_c).length_squared()) / 2) 

        if rmsd > 0.5: continue

        print rmsd

        pose.dump_pdb('debug/test_{0}.pdb'.format(index)) ###DEBUG
        
    exit()
    pose.fold_tree(old_ft)
