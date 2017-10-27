import datetime

import numpy as np

import pyrosetta
from pyrosetta import rosetta

import site_settings
import geometry
import pose_analysis
import assemble_helpers
import IO


def xyz_to_np_array(xyz):
    '''Convert an xyz vector to a numpy array.'''
    return np.array([xyz.x, xyz.y, xyz.z])

def extract_fragment(lh_library, loop_size, frag_index):
    '''Get the fragment from a loophash library given its size and index.'''
    hashmap = lh_library.gethash(loop_size)
   
    backbone_seg = rosetta.protocols.loophash.BackboneSegment()
    bb_data = rosetta.protocols.loophash.BBData()
    extra_data = rosetta.protocols.loophash.BBExtraData()

    # Get the backbone

    cp = hashmap.get_peptide( frag_index )
    lh_library.backbone_database().get_backbone_segment( cp.index, cp.offset , loop_size, backbone_seg )
   
    # Get the sequence

    lh_library.backbone_database().get_protein( cp.index, bb_data )
    lh_library.backbone_database().get_extra_data(bb_data.extra_key, extra_data)

    assert(cp.offset % 3 == 0)
    seq_offset = int(cp.offset/3)
    sequence = extra_data.sequence
    loop_sequence = sequence[seq_offset: seq_offset + loop_size]

    return backbone_seg, loop_sequence

def generate_preproteins(pose, loop, lh_library, ref_seqposes, ref_ca_coords):
    '''Generate preprotein structures given a loop defined as the pair (start, stop).
    The torsions of start and stop will not be changed
    '''
    # Mutate the pose to all ALA
   
    mutater = rosetta.protocols.simple_moves.MutateResidue()
    mutater.set_res_name("ALA")
    for i in range(1, pose.size() + 1):
        mutater.set_target(i)
        mutater.apply(pose)
    
    for i in [33, 34, 43, 44]: ###DEBUG
        pose.real_to_virtual(i)
    
    # Get residues in the two chains

    chain1 = list(range(1, loop[0] + 1))
    chain2 = [i for i in range(loop[1], pose.size() + 1) if not pose.residue(i).is_virtual_residue()]

    # Get the indices of fragments

    loop_size = loop[1] - loop[0] - 1
    hashmap = lh_library.gethash(loop_size)
    
    leap_index_bucket = rosetta.std.vector_unsigned_long()
    #for i in range(10000):
    for i in range(hashmap.n_loops()):
        leap_index_bucket.append(i)

    # Insert the fragments

    for leap_index in leap_index_bucket:
        bb_seg, sequence = extract_fragment(lh_library, loop_size, leap_index)

        phis = bb_seg.phi()
        psis = bb_seg.psi()
        omegas = bb_seg.omega()

        for j in range(len(phis)):
                pose.set_phi(loop[0] + j + 1, phis[j])
                pose.set_psi(loop[0] + j + 1, psis[j])
                pose.set_omega(loop[0] + j + 1, omegas[j])

        ca_coords = [xyz_to_np_array(pose.residue(i).xyz("CA")) for i in ref_seqposes]
        rmsd = geometry.RMSD(ca_coords, ref_ca_coords)
        
        # Filter the results
        
        if rmsd > 10: continue
            
        buhs = pose_analysis.get_buried_unsatisfied_hbonds(pose)
        if sum(buhs[i] for i in range(loop[0] - 1, loop[1] + 2)) > 0 : continue
        
        if not pose_analysis.check_clashes(pose): continue

        inter_chain_distance1 = pose_analysis.calc_ca_distances_between_residue_groups(pose, chain1, chain2)
        inter_chain_distance2 = pose_analysis.calc_ca_distances_between_residue_groups(pose, chain2, chain1)

        if np.mean(inter_chain_distance1) < 9 or np.mean(inter_chain_distance1) > 13 or max(inter_chain_distance1) > 18: continue
        if np.mean(inter_chain_distance2) < 9 or np.mean(inter_chain_distance2) > 13 or max(inter_chain_distance2) > 18: continue

        print leap_index, rmsd, np.mean(inter_chain_distance1), np.mean(inter_chain_distance2), max(inter_chain_distance1), max(inter_chain_distance2)
        
        for i in [33, 34, 43, 44]: ###DEBUG
            pose.virtual_to_real(i)
        pose.dump_pdb('data/preproteins_3_8_helix_20/preprotein_{0}.pdb'.format(leap_index)) #DEBUG
        for i in [33, 34, 43, 44]: ###DEBUG
            pose.real_to_virtual(i)

def append_alanines(pose, num):
    '''Append num alanines to the end of a pose'''
    residue_type_set = pose.residue_type_set_for_pose()
    
    for i in range(num):
        new_rsd = rosetta.core.conformation.ResidueFactory.create_residue( residue_type_set.name_map("ALA") )
        rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(pose, pose.size())
        pose.append_residue_by_bond(new_rsd, True)
        pose.set_omega(pose.size() - 1, 180);

def get_ordered_peptides(pose, connections_seqpos):
    '''Get the start and stop of peptides in the order of the assembled protein.
    Return a list of tuples (start, stop, connection_id) where start and stop points
    to the position of residues in the current pose.
    '''
    pps = []

    # Find the first peptide in the pose

    for i, c in enumerate(connections_seqpos):
        start = pose.fold_tree().boundary_left(c[0])

        no_prev = True
        for cc in connections_seqpos:
            if cc[1] == start:
                no_prev = False

        if no_prev:
            pps.append((start, c[0], i))
  
    # Find subsequent peptides

    while True:
        start = connections_seqpos[pps[-1][2]][1]
        stop = pose.fold_tree().boundary_right(start + 1) ###NOTE: here assumes that there is no 1 residue peptide

        c_id = -1
        for i, c in enumerate(connections_seqpos):
            if c[0] == stop:
                c_id = i

        pps.append((start, stop, c_id))

        if c_id == -1: break

    return pps

def get_connections_after_adding_linkers(connections_seqpos, pps):
    '''Calcuate the connections after adding the linkers given the
    original connections and the ordered peptides.
    '''
    connections_after_adding_linkers = connections_seqpos[:]
    current_position = 1

    for pp in pps:
        if pp[2] == -1: break
        
        current_position += pp[1] - pp[0]
        stop = current_position + connections_seqpos[pp[2]][2] + 1

        connections_after_adding_linkers[pp[2]] = (current_position, stop, connections_seqpos[pp[2]][2])
        current_position = stop

    return connections_after_adding_linkers

def add_linkers(pose, connections, seqpos_map):
    '''Add linkers to the pose and change the order of the pose 
    to the order of the aseembled structure.
    connections is a list of tuples (res1, res2, connection_length) where residues are 
    specified by PDB indexing as (pose_id, chain, res_id). The PDB indices 
    are converted to Rosetta indices by the seqpos_map.
    '''
    connections_seqpos = [(seqpos_map[c[0]], seqpos_map[c[1]], c[2]) for c in connections]
    
    # Get the start and stop of peptides in the order of the assembled protein

    pps = get_ordered_peptides(pose, connections_seqpos)

    # Get connections after linkers are addded

    connections_after_adding_linkers = get_connections_after_adding_linkers(connections_seqpos, pps)

    new_pose = rosetta.core.pose.Pose()

    for i, pp in enumerate(pps):
        seqposes = rosetta.utility.vector1_unsigned_long()
        for seqpos in range(pp[0], pp[1] + 1):
            seqposes.append(seqpos)
        
        pp_pose = rosetta.core.pose.Pose()
        rosetta.core.pose.pdbslice(pp_pose, pose, seqposes)

        if -1 != pp[2]:
            append_alanines(pp_pose, connections_seqpos[pp[2]][2])

        if i == 0:
            new_pose = pp_pose

        else:
            gap_res = new_pose.size()
            new_pose.append_pose_by_jump(pp_pose, new_pose.size())
            rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(new_pose, gap_res + 1)
            new_pose.conformation().declare_chemical_bond(gap_res, "C", gap_res + 1, "N")

    new_pose.pdb_info(rosetta.core.pose.PDBInfo(new_pose.size()))
    new_pose.conformation().reset_chain_endings();

    return new_pose

def pilot_connect_chains():
    pdb_files = ['data/antiparallel_sheets_3_8/2_2_30_30/sheet.pdb',
                 'data/straight_helices/20/helix.pdb']
    transformation_files = ['database/transformations/sheet_helix_transformation.json']
    transformation_residue_pairs = [((0, 'B', 13), (1, 'A', 13))]
    movable_jumps = [3] 
    connections = [((0, 'B', 16), (0, 'A', 1), 2),
                   ((0, 'C', 24), (0,'B', 9), 2),
                   ((1, 'A', 20), (0,'C', 17), 4)]

    poses = []
    for pdb_file in pdb_files:
        poses.append( rosetta.core.pose.Pose() )
        rosetta.core.import_pose.pose_from_file(poses[-1], pdb_file)

    seqpos_map = assemble_helpers.make_seqpos_map(poses)

    # Do the pre-assembly

    Ts = [IO.load_rigid_body_transformation_from_file(tf) for tf in transformation_files]

    assemble_helpers.pre_assemble(poses, transformation_residue_pairs, Ts)
  
    merged_pose = poses[0]

    merged_pose = add_linkers(merged_pose, connections, seqpos_map)

    return merged_pose


def pilot_make_preproteins():
    ss = site_settings.load_site_settings()
    rosetta.basic.options.set_path_option('lh:db_path', ss['loophash_database_path'])
    loopsizes = rosetta.utility.vector1_unsigned_long()
    loopsizes.append(4)

    lh_library = rosetta.protocols.loophash.LoopHashLibrary(loopsizes)
    lh_library.load_mergeddb()

    #pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(pose, 'data/reassemble_test_assemble_3_8_helix_20/0/assembled.pdb')
    
    pose = pilot_connect_chains()
    ref_helix_cas = [xyz_to_np_array(pose.residue(i).xyz("CA")) for i in range(1, 20)]
    
    assemble_helpers.fast_loop_build(pose, [(32, 35), (42, 45)])
    #pose.dump_pdb('debug/test.pdb')###DEBUG
    
    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(pose.size(), 1, -1)
    pose.fold_tree(ft)
    
    rosetta.core.conformation.idealize_position(24, pose.conformation())
   

    loop_start = 20
    loop_stop = loop_start + 5

    start_time = datetime.datetime.now()
   
    generate_preproteins(pose, (loop_start, loop_stop), lh_library,  list(range(1, 20)), ref_helix_cas)

    end_time = datetime.datetime.now()
    run_time = end_time - start_time
    print run_time.total_seconds()


def assemble(pose, unmovable_connections):
    '''Assembe the structure given the pose and unmovable connections.
    Note that the movable connection should be sampled in the preprotein
    generating step.
    '''
    # Do fast design
    
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
    </TASKOPERATIONS>
    <MOVERS>
        <FastDesign name="fastdes" task_operations="limitchi2,ex12,layer_all" clear_designable_residues="1" repeats="5" ramp_down_constraints="0"/>
    </MOVERS>
    ''')
    fast_design = xmlobj.get_mover('fastdes')

    fast_design.apply(pose)

    # Do loop design for the unmovable loops
    
    ss = site_settings.load_site_settings()

    rosetta.basic.options.set_path_option('lh:db_path', ss['loophash_database_path'])
    loopsizes = rosetta.utility.vector1_int()
    for c in unmovable_connections:
        loopsizes.append(c[1] - c[0] + 3)
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
    loops = assemble_helpers.loops_to_rosetta_loops(unmovable_connections)
    
    tf = loop_modeler.get_task_factory()
    restrict_to_loops = rosetta.protocols.toolbox.task_operations.RestrictToLoopsAndNeighbors()
    restrict_to_loops.set_include_neighbors(True)
    restrict_to_loops.set_cutoff_distance(6)
    restrict_to_loops.set_design_neighbors(True)
    restrict_to_loops.set_loops(loops)
    tf.push_back(restrict_to_loops)

    loop_modeler.fullatom_stage().set_temp_cycles(50, False)
    loop_modeler.set_loops(loops)

    loop_modeler.apply(pose)

if __name__ == '__main__':
    pyrosetta.init()

    pilot_make_preproteins()
