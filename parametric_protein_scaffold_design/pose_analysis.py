import json

import pyrosetta
import pyrosetta.rosetta as rosetta


def get_buried_unsatisfied_hbonds(pose):
    '''Get buried unsatisfied hydrogen bonds.
    Return a vector1 of the number of buried unsatisfied Hbonds at each residue. 
    '''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
    pose.update_residue_neighbors()

    calculator = rosetta.protocols.toolbox.pose_metric_calculators.BuriedUnsatisfiedPolarsCalculator('default', 'default')
    buhs = json.loads(calculator.get("residue_bur_unsat_polars", pose))
    
    buhs1 = rosetta.utility.vector1_unsigned_long()
    for buh in buhs:
        buhs1.append(buh)

    return buhs1

def check_clashes(pose):
    '''Return true if there are no clashes in the pose.'''
    scale_factor = 0.36

    for i in range(1, pose.size() + 1):
        res_i = pose.residue(i)
        if res_i.is_virtual_residue(): continue
        
        for j in range(i + 1, pose.size() + 1):
            if i == j or i == j + 1 or i == j - 1:
                continue

            res_j = pose.residue(j)
            if res_j.is_virtual_residue(): continue

            for ai in range(1, res_i.nheavyatoms() + 1):
                for aj in range(1, res_j.nheavyatoms() + 1):

                    vi = res_i.xyz(ai)
                    vj = res_j.xyz(aj)
                    ri = res_i.atom_type(ai).lj_radius()
                    rj = res_j.atom_type(aj).lj_radius()

                    if (vi - vj).length_squared() < scale_factor * ((ri + rj) ** 2):
                        return False

    return True

def calc_ca_distances_between_residue_groups(pose, group1, group2):
    '''Calcualte the Ca distance between groups of residues.
    Return lists of inter group distances of each residues in the first group.
    '''
    distances = []
    for i in group1:
        dmin = float('inf')
        for j in group2:
            d = (pose.residue(i).xyz('CA') - pose.residue(j).xyz('CA')).length()

            if d < dmin:
                dmin = d
        distances.append(dmin)
    return distances

if __name__ == '__main__':
    import sys
   
    pyrosetta.init()

    pdb = sys.argv[1]
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb)

    buhs = get_buried_unsatisfied_hbonds(pose)
    print buhs
