#!/usr/bin/env python2.7

import pyrosetta
import pyrosetta.rosetta as rosetta
import parametric_protein_scaffold_design as PPSD


def ss_prediction_filter(pose, filter_dict):
    '''Filter by the secondary structure prediction.'''
    
    sspf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<SSPrediction name="sspred" use_svm="false" cmd="./dependencies/dependencies/psipred/runpsipred_single" use_probability="false" confidence="1" threshold="0.75"/>')

    filter_dict['SSPrediction']  = sspf.report_sm(pose)
    return sspf.apply(pose)

def pack_stat_filter(pose, filter_dict):
    '''Filter by the packing statistics filter.'''

    psf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<PackStat name="pack" threshold="0.6" confidence="1"/>')
    
    filter_dict['PackStat'] = psf.report_sm(pose)

    # Because the result of the PackStatFilter is stochastic, do not call the apply() function again

    return filter_dict['PackStat'] > 0.6

def holes_filter(pose, filter_dict):
    '''Filter by the holes filter.'''
    rosetta.basic.options.set_file_option('holes:dalphaball', './dependencies/dependencies/DAlpahBall/DAlphaBall.gcc')

    hf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<Holes name="holes" threshold="2.0" confidence="1"/>')

    filter_dict['Holes'] = hf.report_sm(pose)
    return hf.apply(pose)

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


if __name__ == '__main__':
    pyrosetta.init()

    #pdb_file = 'data/test_assemble/0/assembled.pdb' 
    pdb_file = '/home/xingjie/DataBases/PPIAffinityDatabase/Affinity_Benchmark/1EMV_r_u.pdb' 
    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    filter_dict = {}

    #ss_prediction_filter(pose, filter_dict)
    #pack_stat_filter(pose, filter_dict)
    #holes_filter(pose, filter_dict)
    #helix_shape_complementarity_filter(pose, filter_dict)

    print filter_dict
