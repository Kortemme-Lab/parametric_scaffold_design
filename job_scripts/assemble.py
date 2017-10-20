#!/usr/bin/env python2.7
import sys
import os
import json
import datetime

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

from rosetta.protocols.assembly.secondary_structure_assembly import *


def assemble_antiparallel_sheet_3_8_helix_15(data_path, num_jobs, job_id):
    pdb_files = ['data/antiparallel_sheets_3_8/2_2_30_30/sheet.pdb',
                 'data/straight_helices/15/helix.pdb']
    transformation_files = ['database/transformations/sheet_helix_transformation.json']
    transformation_residue_pairs = [((0, 'B', 13), (1, 'A', 7))]
    movable_jumps = [3] 
    connections = [((0, 'B', 16), (0, 'A', 1), 2),
                   ((0, 'C', 24), (0,'B', 9), 2),
                   ((1, 'A', 15), (0,'C', 17), 4)]

    PPSD.assembly.run_batch_assembly(1, data_path, num_jobs, job_id, 
            pdb_files, transformation_files, transformation_residue_pairs, movable_jumps, connections)

if __name__ == '__main__':

    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    pyrosetta.init()

    #PPSD.assembly.save_residue_transformation('data/test_assemble/after_fast_design.pdb', ('C', 34), ('A', 11), 'data/test_assemble/T_C34_A11.json')
   
    pilot_run(data_path, num_jobs, job_id)
