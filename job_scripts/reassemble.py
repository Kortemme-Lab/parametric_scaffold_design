#!/usr/bin/env python2.7
import sys
import os
import json
import datetime

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD


def reassemble_with_linker_fragments(pose, good_fragments, loops_at_movable_jumps, loops_at_unmovable_jumps):
    '''Reassmble a pose given good fragment for linker loops.
    Loops are defined as pairs of (start, stop).
    '''
    pass

def reassemble_with_linker_fragments_from_file(output_path, input_pdb, good_fragments_file, loops_at_movable_jumps, loops_at_unmovable_jumps):
    '''Reassmble a pose given input files.
    Loops are defined as pairs of (start, stop).
    '''
    task_info = {}
    start_time = datetime.datetime.now()

    # Load fragments

    with open(good_fragments_file, 'r') as gf:
        good_fragments = json.load(gf)

    # Do assembly

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)

    reassemble_with_linker_fragments(pose, good_fragments, loops_at_movable_jumps, loops_at_unmovable_jumps)

    pose.dump_pdb(os.path.join(output_path, 'assembled.pdb'))
    PPSD.io.sequence_to_fasta_file(os.path.join(output_path, 'assembled.fasta'), 'assembled', pose.sequence())

    end_time = datetime.datetime.now()
    run_time = end_time - start_time

    # Save the task info

    task_info['sequence'] = pose.sequence()
    task_info['score'] = pose.energies().total_energy()
    task_info['run_time'] = run_time.total_seconds()

    PPSD.io.save_task_info_file(output_path, task_info)

def run_tasks(task_list, num_jobs, job_id):
    '''Run tasks that belongs to the current job thread'''
    for i, t in enumerate(task_list):
        if job_id == i % num_jobs:
            if not os.path.exists(t['output_path']): 
                os.makedirs(t['output_path'])

            reassemble_with_linker_fragments_from_file(t['output_path'], t['input_pdb'], t['good_fragments_file'], 
                    t['loops_at_movable_jumps'], t['loops_at_unmovable_jumps'])

def pilot_run(data_path, num_jobs, job_id):
    output_path = data_path

    task_list = []

    for i in range(1):
        task_list.append({'output_path' : os.path.join(data_path, str(i)),
                          'input_pdb': 'data/test_assemble_3_8_loophashkic_noRepackDesign/697/assembled.pdb',
                          'good_fragments_file': 'data/test_assemble_3_8_loophashkic_noRepackDesign/selected_fragments/fragments.json',
                          'loops_at_movable_jumps':[(15, 20)],
                          'loops_at_unmovable_jumps':[(27,30), (37,40)]})
    
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
