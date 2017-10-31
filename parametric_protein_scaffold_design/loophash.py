import pyrosetta
from pyrosetta import rosetta

import site_settings


def load_lh_library(loop_sizes):
    '''Load the loop hash library for a list of loop sizes.'''
    
    ss = site_settings.load_site_settings()
    rosetta.basic.options.set_path_option('lh:db_path', ss['loophash_database_path'])
    loopsizes_v1 = rosetta.utility.vector1_unsigned_long()
    
    for i in loop_sizes:
        loopsizes_v1.append(i)

    lh_library = rosetta.protocols.loophash.LoopHashLibrary(loopsizes_v1)
    lh_library.load_mergeddb()

    return lh_library

def extract_lh_fragment(lh_library, loop_size, frag_index):
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

def get_linkers_from_loophash_db(lh_library, loop_size):
    '''Get the linker fragments from the loop hash database.'''
    linkers = []
    
    # Get the indices of fragments

    hashmap = lh_library.gethash(loop_size)
    
    leap_index_bucket = rosetta.std.vector_unsigned_long()
    #for i in range(1000):
    for i in range(hashmap.n_loops()):
        leap_index_bucket.append(i)

    # Insert the fragments

    for leap_index in leap_index_bucket:
        bb_seg, sequence = extract_lh_fragment(lh_library, loop_size, leap_index)

        linkers.append({'phis': bb_seg.phi(),
                        'psis': bb_seg.psi(),
                        'omegas' : bb_seg.omega()})
    
    return linkers


