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


if __name__ == '__main__':
    import sys
   
    pyrosetta.init()

    pdb = sys.argv[1]
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb)

    buhs = get_buried_unsatisfied_hbonds(pose)
    print buhs
