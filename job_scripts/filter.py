#!/usr/bin/env python2.7
import os
import json

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD


def ss_prediction_filter(pose, filter_dict):
    '''Filter by the secondary structure prediction.'''
    
    sspf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<SSPrediction name="sspred" use_svm="false" cmd="./dependencies/dependencies/psipred/runpsipred_single" use_probability="false" confidence="1" threshold="0.75"/>')

    filter_dict['SSPrediction']  = {'score' : sspf.report_sm(pose),
                                    'pass' : sspf.apply(pose)}
    return filter_dict['SSPrediction']['pass']

def pack_stat_filter(pose, filter_dict):
    '''Filter by the packing statistics filter.'''

    psf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<PackStat name="pack" threshold="0.6" confidence="1"/>')
    
    score = psf.report_sm(pose)
    filter_dict['PackStat'] = {'score' : score,
                               'pass' : score > 0.6} 

    # Because the result of the PackStatFilter is stochastic, do not call the apply() function again

    return filter_dict['PackStat']['pass']

def holes_filter(pose, filter_dict):
    '''Filter by the holes filter.'''
    rosetta.basic.options.set_file_option('holes:dalphaball', './dependencies/dependencies/DAlpahBall/DAlphaBall.gcc')

    hf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<Holes name="holes" threshold="2.0" confidence="1"/>')

    filter_dict['Holes'] = {'score' : hf.report_sm(pose),
                            'pass' : hf.apply(pose)}
    return filter_dict['Holes']['pass']

def helix_shape_complementarity_filter(pose, filter_dict):
    '''Filter by helix shape complementarity.'''
    #hscf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<SSShapeComplementarity name="sc_hx" />')
    #hscf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<SSShapeComplementarity name="sc_hx" helices="true" loops="false" confidence="1" verbose="1" threshold="0.6" />')

    #filter_dict['SSShapeComplementarity'] = hscf.report_sm(pose)
    #return hscf.apply(pose)

    #ff = rosetta.protocols.filters.FilterFactory.get_instance()
    #fmap = ff.filter_creator_map()
    #
    #for f in fmap:
    #    print f

    print "Not working yet!"

def buried_unsatisfied_hbond_filter(pose, filter_dict):
    '''Filter by the number of buried unsatisfied hydrogen bonds.'''
    buhf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<BuriedUnsatHbonds name="buriedunsat" jump_number="0" cutoff="20" />')

    filter_dict['BuriedUnsatHbonds'] = {'score' : buhf.report_sm(pose),
                                        'pass' : buhf.apply(pose)}
    return filter_dict['BuriedUnsatHbonds']['pass']

def filter_one_design(design_path):
    '''Filter one design inside the design path. And save the filter information
       inside the filter_info.json file.
    '''
    # Load the design

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, os.path.join(design_path, 'assembled.pdb'))

    # Filter the design

    filter_dict = {}

    filter_list = [ss_prediction_filter, pack_stat_filter, 
            holes_filter, buried_unsatisfied_hbond_filter]

    for f in filter_list:
        f(pose, filter_dict)

    # Dump the results

    with open(os.path.join(design_path, 'filter_info.json'), 'w') as f:
        json.dump(filter_dict, f)

def select_designs(input_path, output_path, max_pass):
    '''Select the designs in the input path and link the designs that pass the filter to the output_path.
    If the number of passed designs are more than the max_pass, only save the top max_pass designs.
    '''
    # Load the designs

    designs = []

    for d in os.listdir(input_path):
        try:
            with open(os.path.join(input_path, d, 'task_info.json'), 'r') as f:
                task_info = json.load(f)

            with open(os.path.join(input_path, d, 'filter_info.json'), 'r') as f:
                filter_info = json.load(f)

        except IOError:
            continue

        designs.append({'id':d,
                        'path':os.path.join(input_path, d),
                        'task_info':task_info,
                        'filter_info':filter_info})
   
    designs = sorted(designs, key=lambda x : x['task_info']['score'])

    #TODO

if __name__ == '__main__':
    pyrosetta.init()

    #pdb_file = 'data/test_assemble/0/assembled.pdb' 
    #pdb_file = '/home/xingjie/DataBases/PPIAffinityDatabase/Affinity_Benchmark/1EMV_r_u.pdb' 
    #pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    #filter_dict = {}

    #ss_prediction_filter(pose, filter_dict)
    #pack_stat_filter(pose, filter_dict)
    #holes_filter(pose, filter_dict)
    #helix_shape_complementarity_filter(pose, filter_dict)
    #buried_unsatisfied_hbond_filter(pose, filter_dict)

    #print filter_dict

    #filter_designs('data/test_assemble', 'data/test_filter', 1)

    filter_one_design('data/test_assemble/0')
