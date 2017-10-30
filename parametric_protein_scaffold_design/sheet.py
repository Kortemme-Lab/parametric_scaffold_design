import numpy as np
import cylinder_fitting

import pyrosetta
from pyrosetta import rosetta


def fit_sheet_to_cylinder(pose, strands, bonding_pairs):
    '''Fit a beta sheet to a cylinder. The sheet is defined by a 
    list of strands. Each strand is a list of seqposes. The difference
    between two consecutive seqposes in a strand should be 2. For example
    [5, 7, 9].
    During the fitting, extra points will be added to consecutive residues
    and at the middle of bonding pairs to avoid bad fitting. The bonding pairs
    are pairs of ((strand_id, residue_id), (strand_id, residue_id)).
    Return the radius of the cylinder and the crease strand angle in radians
    at the middle seqpos of the middle strand.
    '''
    points = []

    for strand in strands:
        for i, seqpos in enumerate(strand):
            xyz = pose.residue(seqpos).xyz('CA')
            p = np.array([xyz.x, xyz.y, xyz.z])
            points.append(p)

            if i < len(strand) - 1:
                xyz_next = pose.residue(strand[i + 1]).xyz('CA')
                p_next = np.array([xyz_next.x, xyz_next.y, xyz_next.z])  
                points.append((p + p_next) / 2)

    for r1, r2 in bonding_pairs:
        xyz1 = pose.residue(strands[r1[0]][r1[1]]).xyz('CA')
        p1 = np.array([xyz1.x, xyz1.y, xyz1.z])
        xyz2 = pose.residue(strands[r2[0]][r2[1]]).xyz('CA')
        p2 = np.array([xyz2.x, xyz2.y, xyz2.z])

        points.append((p1 + p2) / 2)

    # Fit the sheet into a cylinder
    
    direction, center, radius, fit_err = cylinder_fitting.fit(points)

    # Get the strand direction at the center of the central strand

    cs = int(len(strands) / 2)
    cc = int(len(strands[cs]) / 2)

    cs_dir_xyz = pose.residue(strands[cs][cc + 1]).xyz('CA') - pose.residue(strands[cs][cc]).xyz('CA')
    cs_dir = np.array([cs_dir_xyz.x, cs_dir_xyz.y, cs_dir_xyz.z])

    # Calculate the angle between the strand direction and cylinder direction

    cs_p_xyz = pose.residue(strands[cs][cc]).xyz('CA')
    cs_p = np.array([cs_p_xyz.x, cs_p_xyz.y, cs_p_xyz.z]) - center

    sin_abs = np.linalg.norm(np.cross(cs_dir, direction)) / (np.linalg.norm(cs_dir) * np.linalg.norm(direction))
    sign = np.sign(np.dot(cs_dir, direction) * np.dot(cs_dir, np.cross(direction, cs_p)))
    crease_strand_angle = sign * np.arcsin(sin_abs)

    return radius, crease_strand_angle

def fit_sheet_to_cylinder_pdb(pdb_file, strands_pdb, bonding_pairs=[]):
    '''Fit a beta sheet from a PDB file to a cylinder. The sheet id defined
    by a list of strands. Each strand is a list of pbd index pairs. For example
    [('A', 3), ('A', 5), ('A', 7)].
    During the fitting, extra points will be added to consecutive residues
    and at the middle of bonding pairs to avoid bad fitting. The bonding pairs
    are pairs of (strand_id, residue_id).
    Return the radius of the cylinder and the crease strand angle in radians
    at the middle seqpos of the middle strand.
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    
    strands = [[pose.pdb_info().pdb2pose(i[0], i[1]) for i in s_pdb ] for s_pdb in strands_pdb]

    return fit_sheet_to_cylinder(pose, strands, bonding_pairs)


if __name__ == '__main__':
    pyrosetta.init()

    for i in range(5):
        for j in range(5):
            d = '{0}_{1}_{2}0_{3}0'.format(i, j, i + 1, j+ 1)

            r, theta = fit_sheet_to_cylinder_pdb('data/antiparallel_sheets_3_8/{0}/sheet.pdb'.format(d), 
                    [[('A', 1), ('A', 3), ('A', 5), ('A', 7)],
                     [('B', 10), ('B', 12), ('B', 14), ('B', 16)],
                     [('C', 17), ('C', 19), ('C', 21), ('C', 23)]])

            print d, r, np.degrees(theta)
