#!/usr/bin/env python2.7
import sys
import os
import json
import datetime

import numpy as np

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD

def assemble_from_file(output_path, preprotein_path, loops):
    '''Reassmble a pose given input files.
    Loops are defined as pairs of (start, stop).
    '''
    task_info = {}
    start_time = datetime.datetime.now()

    # Choose a random preprotein
    
    pp_pdbs = [f for f in os.listdir(preprotein_path) if f.endswith('.pdb')]
    preprotein_pdb = os.path.join(preprotein_path, pp_pdbs[np.random.randint(0, len(pp_pdbs))]) 

    # Do assembly

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, preprotein_pdb)

    PPSD.SSAssembler2.assemble(pose, loops)

    pose.dump_pdb(os.path.join(output_path, 'assembled.pdb'))
    PPSD.IO.sequence_to_fasta_file(os.path.join(output_path, 'assembled.fasta'), 'assembled', pose.sequence())

    end_time = datetime.datetime.now()
    run_time = end_time - start_time

    # Save the task info

    task_info['sequence'] = pose.sequence()
    task_info['score'] = pose.energies().total_energy()
    task_info['run_time'] = run_time.total_seconds()

    PPSD.IO.save_task_info_file(output_path, task_info)

def run_tasks(task_list, num_jobs, job_id):
    '''Run tasks that belongs to the current job thread'''
    for i, t in enumerate(task_list):
        if job_id == i % num_jobs:
            if not os.path.exists(t['output_path']): 
                os.makedirs(t['output_path'])

            assemble_from_file(t['output_path'], t['preprotein_path'], t['loops'])

def pilot_run(data_path, num_jobs, job_id):
    output_path = data_path

    task_list = []

    for i in range(1000):
        task_list.append({'output_path' : os.path.join(data_path, str(i)),
                          #'preprotein_path': 'data/preproteins_3_8_helix_20',
                          'preprotein_path': 'data/preproteins_3_8_helix_20_new',
                          'loops':[(32,35), (42,45)]})
    
    run_tasks(task_list, num_jobs, job_id)


if __name__ == '__main__':
    
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    pyrosetta.init()
    
    pilot_run(data_path, num_jobs, job_id)
