#!/usr/bin/env python2.7
import sys
import os
import json
import datetime

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD


def refine(pose, repack_positions=[]):
    '''Refine a design with fast design.
    Use the repack_positions argument to restrict positions that are only repackable.
    '''
   
    # Make the FastDesign mover
    restrict_repacking_operation = \
    '''
    <OperateOnResidueSubset name="restrict_repacking">
        <Index resnums="{0}" />
        <RestrictToRepackingRLT />
    </OperateOnResidueSubset>'''.format(','.join(str(p) for p in repack_positions)) if len(repack_positions) > 0 else ''


    xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <TASKOPERATIONS>
        <LimitAromaChi2 name="limitchi2" include_trp="1" />
        <ExtraRotamersGeneric name="ex12" ex1="true" ex2="true" />
        <LayerDesign name="layer_all" layer="core_boundary_surface_Nterm_Cterm" use_sidechain_neighbors="True" pore_radius="2.0" verbose="true" core="3.5" surface="0.95" >
    		<Nterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Nterm>
    		<Cterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Cterm>
        </LayerDesign>
        {0}
    </TASKOPERATIONS>
    <MOVERS>
        <FastDesign name="fastdes" task_operations="limitchi2,ex12,layer_all{1}" clear_designable_residues="1" repeats="5" ramp_down_constraints="0" />
    </MOVERS>
    '''.format(restrict_repacking_operation, ',restrict_repacking' if len(repack_positions) > 0 else ''))
    fast_design = xmlobj.get_mover('fastdes')

    # Refine
    
    fast_design.apply(pose)

def refine_from_file(input_pdb, repack_positions, output_path):
    '''Refine a design given the input_pdb file and the positions restricted for repacking'''
    task_info = {}
    start_time = datetime.datetime.now()
    
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, input_pdb)

    refine(pose, repack_positions)

    pose.dump_pdb(os.path.join(output_path, 'assembled.pdb'))
    PPSD.IO.sequence_to_fasta_file(os.path.join(output_path, 'assembled.fasta'), 'assembled', pose.sequence())
    
    end_time = datetime.datetime.now()
    run_time = end_time - start_time
    
    # Save task information
    
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

            refine_from_file(t['input_pdb'], t['repack_positions'], t['output_path'])

def pilot_run(data_path, num_jobs, job_id):
    input_pdb = 'data/test_assemble_3_8_restrict_turns/961/assembled.pdb'
    output_path = data_path
    task_list = []

    for i in range(1000):
        task_list.append({'input_pdb':input_pdb,
                          'repack_positions' : [],# [15,16,17,18,19,20,27,28,29,30,37,38,39,40],
                          'output_path':os.path.join(output_path, str(i))})
    
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
