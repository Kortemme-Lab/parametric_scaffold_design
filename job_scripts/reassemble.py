#!/usr/bin/env python2.7
import sys
import os
import json
import datetime

import numpy as np

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD


def reassemble_with_linker_fragments(pose, good_fragments, loops):
    '''Reassmble a pose given good fragment for linker loops.
    Loops are defined as pairs of (start, stop).
    '''

    def find_loop_id(start, stop):
        '''Find the id of the loop in the good_fragments'''
        for i in range(len(good_fragments)):
            if good_fragments[i]['start'] == start \
               and good_fragments[i]['stop'] == stop \
               and len(good_fragments[i]['fragments']) > 0:
                return i
       
        raise Exception("Cannot find fragments for loop ({0}, {1})".format(start, stop))

    # Mutate the pose to all ALA
   
    mutater = rosetta.protocols.simple_moves.MutateResidue()
    mutater.set_res_name("ALA")
    for i in range(1, pose.size() + 1):
        mutater.set_target(i)
        mutater.apply(pose)

    # Fragment insertion for loops

    for l in loops:
        loop_id = find_loop_id(l[0], l[1]) 
        frag = good_fragments[loop_id]['fragments'][np.random.randint(0, len(good_fragments[loop_id]['fragments']))]

        #print frag['design_id'] ###DEBUG
        #print frag ###DEBUG

        for i in range(len(frag['sequence'])):
            seqpos = l[0] + i
            
            # Mutate the sequence
            
            mutater.set_res_name(frag['sequence'][i])
            mutater.set_target(seqpos)
            mutater.apply(pose)

            # Change the torsions

            pose.set_phi(seqpos, frag['torsions'][3 * i])
            pose.set_psi(seqpos, frag['torsions'][3 * i + 1])
            pose.set_omega(seqpos, frag['torsions'][3 * i + 2])

    # Do fast design while keep the inserted fragments fixed    

    loop_positions = [str(i) for l in loops for i in range(l[0], l[1] + 1)]

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
        <OperateOnResidueSubset name="restrict_repacking">
            <Index resnums="{0}" />
            <RestrictToRepackingRLT />
        </OperateOnResidueSubset> 
    </TASKOPERATIONS>
    <MOVERS>
        <FastDesign name="fastdes" task_operations="limitchi2,ex12,layer_all,restrict_repacking" clear_designable_residues="1" repeats="5" ramp_down_constraints="0"/>
    </MOVERS>
    '''.format(','.join(loop_positions)))
    fast_design = xmlobj.get_mover('fastdes')

    fast_design.apply(pose)

def reassemble_with_linker_fragments_from_file(output_path, input_pdb, good_fragments_file, loops):
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

    reassemble_with_linker_fragments(pose, good_fragments, loops)

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

            reassemble_with_linker_fragments_from_file(t['output_path'], t['input_pdb'], t['good_fragments_file'], 
                    t['loops'])

def pilot_run(data_path, num_jobs, job_id):
    output_path = data_path

    task_list = []

    for i in range(1000):
        task_list.append({'output_path' : os.path.join(data_path, str(i)),
                          #'input_pdb': 'data/test_assemble_3_8_loophashkic_noBBRelax/0/assembled.pdb',
                          #'good_fragments_file': 'data/test_assemble_3_8_loophashkic_noBBRelax/selected_fragments/fragments.json',
                          #'loops':[(15, 20), (27,30), (37,40)]})
                          #'input_pdb': 'data/test_assemble_3_8_helix_20/0/assembled.pdb',
                          #'good_fragments_file': 'data/test_assemble_3_8_helix_20/selected_fragments/fragments.json',
                          'input_pdb': 'data/test_assemble_3_8_helix_20_LeuDock/0/assembled.pdb',
                          'good_fragments_file': 'data/test_assemble_3_8_helix_20_LeuDock/selected_fragments/fragments.json',
                          'loops':[(20, 25), (32,35), (42,45)]})
    
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
