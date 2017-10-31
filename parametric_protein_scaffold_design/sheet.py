import numpy as np
import cylinder_fitting

import pyrosetta
from pyrosetta import rosetta


def fit_sheet_to_cylinder(pose, strands):
    '''Fit a beta sheet to a cylinder. The sheet is defined by a 
    list of strands. Each strand is a list of seqposes. The difference
    between two consecutive seqposes in a strand should be 2. For example
    [5, 7, 9].
    Return the radius of the cylinder and the crease strand angle in radians
    at the middle seqpos of the middle strand.
    '''
    points = []

    for strand in strands:
        for i, seqpos in enumerate(strand):
            for atom in ['CA', 'C', 'N', 'O']:
                xyz = pose.residue(seqpos).xyz(atom)
                p = np.array([xyz.x, xyz.y, xyz.z])
                points.append(p)

    # Fit the sheet into a cylinder
    
    direction, center, radius, fit_err = cylinder_fitting.fit(points)

    #cylinder_fitting.show_fit(direction, center, radius, points) ###DEBUG

    # Get the strand direction at the center of the central strand

    cs = int(len(strands) / 2)
    cc = int(len(strands[cs]) / 2)

    cs_dir_xyz = pose.residue(strands[cs][cc + 2]).xyz('CA') - pose.residue(strands[cs][cc]).xyz('CA')
    cs_dir = np.array([cs_dir_xyz.x, cs_dir_xyz.y, cs_dir_xyz.z])

    # Calculate the angle between the strand direction and cylinder direction

    cs_p_xyz = pose.residue(strands[cs][cc]).xyz('CA')
    cs_p = np.array([cs_p_xyz.x, cs_p_xyz.y, cs_p_xyz.z]) - center

    sin_abs = np.linalg.norm(np.cross(cs_dir, direction)) / (np.linalg.norm(cs_dir) * np.linalg.norm(direction))
    sign = np.sign(np.dot(cs_dir, direction) * np.dot(cs_dir, np.cross(direction, cs_p)))
    crease_strand_angle = sign * np.arcsin(sin_abs)

    return radius, crease_strand_angle

def fit_sheet_to_cylinder_pdb(pdb_file, strands_pdb):
    '''Fit a beta sheet from a PDB file to a cylinder. The sheet id defined
    by a list of strands. Each strand is a tuple of (chain, start, stop)
    During the fitting, extra points will be added to consecutive residues
    and at the middle of bonding pairs to avoid bad fitting. The bonding pairs
    are pairs of (strand_id, residue_id).
    Return the radius of the cylinder and the crease strand angle in radians
    at the middle seqpos of the middle strand.
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    
    strands = [[pose.pdb_info().pdb2pose(s_pdb[0], i) for i in range(s_pdb[1], s_pdb[2] + 1)] for s_pdb in strands_pdb]

    return fit_sheet_to_cylinder(pose, strands)


def calc_square_sheet_bending_angles(pose, center, corner1, corner2, corner3, corner4):
    '''Cacluate the two bending angles of a square sheet.'''
    def get_position(seqpos):
        xyz = pose.residue(seqpos).xyz('CA')
        return np.array([xyz.x, xyz.y, xyz.z])

    def angle(vector1, vector2):
        return np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))

    pc = get_position(center)
    pco1 = get_position(corner1)
    pco2 = get_position(corner2)
    pco3 = get_position(corner3)
    pco4 = get_position(corner4)
    pm12 = get_position(int((corner1 + corner2) / 2))
    pm34 = get_position(int((corner3 + corner4) / 2))

    v1 = pco1 - pc
    v2 = pco2 - pc
    v3 = pco3 - pc
    v4 = pco4 - pc

    a12 = angle(pco1 - pm12, pco2 - pm12)
    a34 = angle(pco3 - pm34, pco4 - pm34)

    return angle(v1, v3), angle(v2, v4), (a12 + a34) / 2

def calc_square_sheet_bending_angles_pdb(pdb_file, center, corner1, corner2, corner3, corner4):
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    c = pose.pdb_info().pdb2pose(center[0], center[1])
    co1 = pose.pdb_info().pdb2pose(corner1[0], corner1[1])
    co2 = pose.pdb_info().pdb2pose(corner2[0], corner2[1])
    co3 = pose.pdb_info().pdb2pose(corner3[0], corner3[1])
    co4 = pose.pdb_info().pdb2pose(corner4[0], corner4[1])

    return calc_square_sheet_bending_angles(pose, c, co1, co2, co3, co4)

if __name__ == '__main__':
    pyrosetta.init()

    for i in range(5):
        for j in range(5):
            d = '{0}_{1}_{2}0_{3}0'.format(i, j, i + 1, j+ 1)

            #r, theta = fit_sheet_to_cylinder_pdb('data/antiparallel_sheets_3_8/{0}/sheet.pdb'.format(d), 
            #            [('A', 1, 8), ('B', 9, 16), ('C', 17, 24)])

            #print d, r, np.degrees(theta)
            a1, a2, a3 = calc_square_sheet_bending_angles_pdb('data/antiparallel_sheets_3_8/{0}/sheet.pdb'.format(d),
                    ('B', 12), ('C', 17), ('C', 24), ('A', 8), ('A', 1))

            print d, np.degrees(a1), np.degrees(a2), np.degrees(a3)

    designs = ['data/assemble2_3_8_helix_20_new/580/assembled.pdb',
               'data/assemble2_3_8_helix_20_new/94/assembled.pdb',
               'data/assemble2_3_8_helix_20_new/253/assembled.pdb',
               'data/assemble2_3_8_helix_20_new/29/assembled.pdb',
               'data/assemble2_3_8_helix_20_new/733/assembled.pdb',
               '/home/xingjie/Softwares/scripts/forward_folding/data/assemble2_3_8_helix_20_new_580/outputs/pdbs/1477_S_00000009.pdb.gz',
               '/home/xingjie/Softwares/scripts/forward_folding/data/assemble2_3_8_helix_20_new_94/outputs/pdbs/193_S_00000003.pdb.gz',
               '/home/xingjie/Softwares/scripts/forward_folding/data/assemble2_3_8_helix_20_new_253/outputs/pdbs/741_S_00000002.pdb.gz',
               '/home/xingjie/Softwares/scripts/forward_folding/data/assemble2_3_8_helix_20_new_29/outputs/pdbs/848_S_00000007.pdb.gz',
               '/home/xingjie/Softwares/scripts/forward_folding/data/assemble2_3_8_helix_20_new_733/outputs/pdbs/1301_S_00000001.pdb.gz']

    for d in designs:
        #r, theta = fit_sheet_to_cylinder_pdb(d,
        #        [('A', 25, 32), ('A', 35, 42), ('A', 45, 52)])
        #print d
        #print r, np.degrees(theta)

        a1, a2, a3 = calc_square_sheet_bending_angles_pdb(d,
                ('A', 38), ('A', 25), ('A', 32), ('A', 52), ('A', 45))

        print np.degrees(a1), np.degrees(a2), np.degrees(a3)
