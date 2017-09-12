#!/usr/bin/env python2.7
import sys
import os

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

from rosetta.protocols.parametric_design import *


def make_straight_helices(data_path, lengthes):
    '''Make a set of straight helices with different lengthes.'''
    
    for length in lengthes:
        helix_path = os.path.join(data_path, str(length))

        if not os.path.exists(helix_path):
            os.mkdir(helix_path)

        make_a_straight_helix(helix_path, length)

def make_a_straight_helix(helix_path, helix_length):
    '''Make a straight helix and save it to a PDB file.'''

    phd = ParametricHelixDesigner()

    directions = pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    
    for i in range(helix_length - 2):
        directions.append(rosetta.numeric.xyzVector_double_t(0, 0, 1))

    phd.set_directions(directions)
    phd.set_phase_vector(rosetta.numeric.xyzVector_double_t(1, 0, 0))
    phd.set_constraint_weight(1)
    
    helix = rosetta.core.pose.Pose()
    phd.apply(helix)

    helix.dump_pdb(os.path.join(helix_path, 'helix.pdb'))

if __name__ == '__main__':
    
    data_path = sys.argv[1]
    
    pyrosetta.init()

    make_straight_helices(data_path, [10, 20]) 

