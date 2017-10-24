import sys
import os
import json
import datetime

import pyrosetta
import pyrosetta.rosetta as rosetta

import IO
import site_settings
from rosetta.protocols.assembly.secondary_structure_assembly import *


def print_pose_pdb_map(pose):
    '''Print the map between pose indexing and the pdb indexing.'''
    for i in range(1, pose.size() + 1):
        print i, pose.pdb_info().chain(i), pose.pdb_info().number(i)
    
    print '\n'

def save_residue_transformation(pdb_file, res1, res2, fout):
    '''Save the rigid body transformation from the res1 frame to the res2 frame.
    Residues are specified by a tuple of (chain_id, residue_id).
    '''
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    
    seqpos1 = pose.pdb_info().pdb2pose(res1[0], res1[1])
    seqpos2 = pose.pdb_info().pdb2pose(res2[0], res2[1])

    T = calc_residue_frame_transformation(pose, seqpos1, seqpos2)
    IO.write_rigid_body_transformation_to_file(T, fout)

def make_seqpos_map(poses):
    '''Make a map from the (structure_id, chain, pdb_number)
    to the merged pose. The structure_id starts from 0.
    '''
    seqpos_map = {}

    infos = [pose.pdb_info() for pose in poses]

    offset = 0
    for i in range(len(poses)):
        for j in range(1, poses[i].size() + 1):
            seqpos_map[(i, str(infos[i].chain(j)), int(infos[i].number(j)))] = j + offset
        offset += poses[i].size()

    return seqpos_map

def pre_assemble(poses, res_pairs, frame_transformations):
    '''Assemble a list of poses into one pose given the relative
    orientation of one residue from each structures. The resulting structure
    will be saved in the first pose.
    
    Residues are specified by a pair of tuples of (pose_id, chain_id, residue_id).
    '''

    def res_to_seqpos(res, poses_before_assemble):
        '''Convert a residue into the seqpos of the assembled structure'''
        offset = 0
        for i in range(res[0]):
            offset += poses_before_assemble[i].size()

        return offset + poses_before_assemble[res[0]].pdb_info().pdb2pose(res[1], res[2])
   
    # Get the seqposes after assembly
    
    seqpos_pairs = []
    
    for res1, res2 in res_pairs:
        seqpos_pairs.append((res_to_seqpos(res1, poses), res_to_seqpos(res2, poses)))

    # Merge poses into pose1

    for i in range(1, len(poses)):
        poses[0].append_pose_by_jump(poses[i], poses[0].size())
   
    # Transform the chains

    for i in range(len(seqpos_pairs)):
        transform_chain(poses[0], seqpos_pairs[i][0], seqpos_pairs[i][1], frame_transformations[i])

def assemble(pose, movable_jumps, connections, seqpos_map, task_info, sasa_threshold=600):
    '''Assemble the secondary structures into a connected protein. 
    movable_jumps is a list of jump ids,
    connections is a list of tuples (res1, res2, connection_length) where residues are 
    specified by PDB indexing as (pose_id, chain, res_id). The PDB indices 
    are converted to Rosetta indices by the seqpos_map.
    '''
    # Turn off the tracers that dump a lot of outputs

    tracers = rosetta.utility.vector1_string()
    tracers.append('protocols.loops.loops_main')
    rosetta.basic.options.set_string_vector_option('out:mute', tracers)
    
    ssa = SecondaryStructureAssembler()

    ssa.set_docking_res_type('LEU')

    ###DEBUG
    #ssa.pass_dock_low_res(True)
    #ssa.pass_dock_high_res(True)
    #ssa.pass_build_loops(True)
    #ssa.pass_fast_design(True)
    #ssa.pass_design_loops(True)
    ####DEBUG 

    for mj in movable_jumps:
        ssa.add_movable_jump(mj)

    for c in connections:
        ssa.add_connection(seqpos_map[c[0]], seqpos_map[c[1]], c[2])

    # Get the connections (loops) after assembly

    connections_after_loop_building = ssa.connections_after_loop_building(pose)

    # Make docking filters
    
    sasas = []
    for mj in movable_jumps:
        sasas.append(rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<Sasa name="sasa" threshold="{0}" jump="{1}" confidence="1"/>'.format(
            sasa_threshold, str(mj))))
   
    for sasa in sasas:
        ssa.add_docking_filter(sasa)

    # Make a customized fast design mover

    xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <TASKOPERATIONS>
        <LimitAromaChi2 name="limitchi2" include_trp="1" />
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
    </TASKOPERATIONS>
    <MOVERS>
        <FastDesign name="fastdes" task_operations="limitchi2,layer_all" clear_designable_residues="0" repeats="5" ramp_down_constraints="0" >
            <MoveMap bb="false" chi="true" jump="false" />
        </FastDesign>
    </MOVERS>
    ''')
    fast_design = xmlobj.get_mover('fastdes')
    ssa.set_fast_design_mover(fast_design)

    # Make a customized loop design mover for DEBUG
    
    ss = site_settings.load_site_settings()

    rosetta.basic.options.set_path_option('lh:db_path', ss['loophash_database_path'])
    loopsizes = rosetta.utility.vector1_int()
    for c in connections_after_loop_building:
        loopsizes.append(c.loop_length + 4)
    rosetta.basic.options.set_integer_vector_option('lh:loopsizes', loopsizes)

    xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <TASKOPERATIONS>
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
    </TASKOPERATIONS>
    <MOVERS>
        <LoopModeler name="loop_modeler" task_operations="layer_all" config="loophash_kic" loophash_perturb_sequence="true" >
            <Build skip="true"/>
            <Centroid skip="true"/>
        </LoopModeler>
    </MOVERS>
    '''
    )
    loop_modeler = xmlobj.get_mover('loop_modeler')
    loops = rosetta.protocols.loops.Loops()
    
    for c in connections_after_loop_building:
        loops.add_loop(c.res1, c.res2)
    
    tf = loop_modeler.get_task_factory()
    restrict_to_loops = rosetta.protocols.toolbox.task_operations.RestrictToLoopsAndNeighbors()
    restrict_to_loops.set_include_neighbors(True)
    restrict_to_loops.set_cutoff_distance(6)
    restrict_to_loops.set_design_neighbors(True)
    restrict_to_loops.set_loops(loops)
    tf.push_back(restrict_to_loops)

    loop_modeler.fullatom_stage().set_temp_cycles(50, False)

    ssa.set_loop_designer(loop_modeler)

    ssa.apply(pose)

    #for i, sasa in enumerate(sasas): #DEBUG
    #    task_info['sasa_{0}'.format(movable_jumps[i])] = sasa.report_sm(pose)

    #pose.energies().show() ###DEBUG

def assemble_from_files(pdb_files, transformation_files, transformation_residue_pairs, movable_jumps, connections, output_path, sasa_threshold=600):
    '''Assemble two secondary structures given the PDB files, transformation
    file, transformation residues, movable jumps and connections.
    The outputs will be saved to the output_path.
    '''
    task_info = {}
    start_time = datetime.datetime.now()

    # Load the structures

    poses = []
    for pdb_file in pdb_files:
        poses.append( rosetta.core.pose.Pose() )
        rosetta.core.import_pose.pose_from_file(poses[-1], pdb_file)

    seqpos_map = make_seqpos_map(poses)

    # Do the pre-assembly

    Ts = [IO.load_rigid_body_transformation_from_file(tf) for tf in transformation_files]

    pre_assemble(poses, transformation_residue_pairs, Ts)
  
    merged_pose = poses[0]

    # Do the assembly

    assemble(merged_pose, movable_jumps, connections, seqpos_map, task_info, sasa_threshold)

    merged_pose.dump_pdb(os.path.join(output_path, 'assembled.pdb'))
    IO.sequence_to_fasta_file(os.path.join(output_path, 'assembled.fasta'), 'assembled', merged_pose.sequence())

    end_time = datetime.datetime.now()
    run_time = end_time - start_time

    # Save the task info

    task_info['sequence'] = merged_pose.sequence()
    task_info['score'] = merged_pose.energies().total_energy()
    task_info['run_time'] = run_time.total_seconds()

    IO.save_task_info_file(output_path, task_info)

def run_tasks(task_list, num_jobs, job_id):
    '''Run tasks that belongs to the current job thread'''
    for i, t in enumerate(task_list):
        if job_id == i % num_jobs:
            if not os.path.exists(t['output_path']): 
                os.makedirs(t['output_path'])

            assemble_from_files(t['pdb_files'], t['transformation_files'], t['transformation_residue_pairs'],
                    t['movable_jumps'], t['connections'], t['output_path'], t['sasa_threshold'])

def run_batch_assembly(num_structures, data_path, num_jobs, job_id, pdb_files, transformation_files, transformation_residue_pairs, 
        movable_jumps, connections, sasa_threshold=600):
    '''Make a batch of assembly'''
    output_path = data_path

    task_list = []

    for i in range(num_structures):
        task_list.append({'pdb_files': pdb_files,
                      'transformation_files': transformation_files,
                      'transformation_residue_pairs':transformation_residue_pairs,
                      'movable_jumps':movable_jumps,
                      'connections':connections,
                      'sasa_threshold':sasa_threshold,
                      'output_path':os.path.join(output_path, str(i))})
    
    run_tasks(task_list, num_jobs, job_id)
    

