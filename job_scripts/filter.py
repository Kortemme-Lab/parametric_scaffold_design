#!/usr/bin/env python2.7
import os
import sys
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
    hscf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<SSShapeComplementarity name="sc_hx" helices="true" loops="false" confidence="1" verbose="1" min_sc="0.6" />')

    filter_dict['SSShapeComplementarity'] = {'score' : hscf.report_sm(pose),
                                             'pass' : hscf.apply(pose)}
    return filter_dict['SSShapeComplementarity']['pass']

def buried_unsatisfied_hbond_filter(pose, filter_dict):
    '''Filter by the number of buried unsatisfied hydrogen bonds.'''
    buhf = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<BuriedUnsatHbonds name="buriedunsat" jump_number="0" cutoff="20" />')

    filter_dict['BuriedUnsatHbonds'] = {'score' : buhf.report_sm(pose),
                                        'pass' : buhf.apply(pose)}
    return filter_dict['BuriedUnsatHbonds']['pass']

def fragment_analysis(design_path):
    '''Do fragment analysis on a design.
    Return the best CRMSD at the worst position.
    '''
    fqa = PPSD.fragment_quality_analysis.FragmentQualityAnalyzer(
            './dependencies/dependencies/psipred/runpsipred_single', 
            #'fragment_picker.linuxclangrelease', 
            '/netapp/home/rpac/Rosetta_57781/main/source/bin/fragment_picker.default.linuxgccrelease',
            #'database/fragment_quality_analysis/small.vall.gz', 
            '/netapp/home/klabqb3backrub/tools/fragment_generation/vall.jul19.2011',
            'database/fragment_quality_analysis/simple.wghts',
            rosetta_database='/netapp/home/rpac/Rosetta_57781/main/database')

    fdf = fqa.pick_fragments(
            os.path.join(design_path, 'assembled.pdb'), 
            os.path.join(design_path, 'assembled.fasta'), 
            design_path)

    crmsds = PPSD.fragment_quality_analysis.FragmentQualityAnalyzer.get_position_crmsd(fdf)

    return max(crmsds)

def filter_one_design(design_path):
    '''Filter one design inside the design path. And save the filter information
       inside the filter_info.json file.
    '''
    # Load the design
    
    if not os.path.exists(os.path.join(design_path, 'assembled.pdb')):
        return

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, os.path.join(design_path, 'assembled.pdb'))

    # Filter the design

    filter_dict = {}

    filter_list = [ss_prediction_filter, pack_stat_filter, 
            holes_filter, helix_shape_complementarity_filter,
            buried_unsatisfied_hbond_filter]

    for f in filter_list:
        f(pose, filter_dict)

    # Do fragment analysis

    worst_crmsd = fragment_analysis(design_path) 
    threshold = 1

    filter_dict['fragment_analysis'] = {'worst_crmsd' : worst_crmsd,
                                        'pass' : worst_crmsd < threshold} 

    # Dump the results

    with open(os.path.join(design_path, 'filter_info.json'), 'w') as f:
        json.dump(filter_dict, f)

def filter_designs(designs_path, num_jobs, job_id):
    '''Filter all designs in the designs_path'''
    design_pathes = [os.path.join(designs_path, p) for p in os.listdir(designs_path)]

    for i, design_path in enumerate(design_pathes):
        if i % num_jobs == job_id:
            filter_one_design(design_path)

def load_designs(input_path):
    '''Load information of designs into a list.'''
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
   
    return sorted(designs, key=lambda x : x['task_info']['score'])

def select_designs(input_path, max_pass):
    '''Select the designs in the input path that pass all the filters.
    If the number of passed designs are more than the max_pass, only save the top max_pass designs.
    '''
    # Load the designs

    designs = load_designs(input_path)
    print 'Number of designs = ', len(designs) 

    # Get designs that pass all the filters  

    passed_designs = []

    for d in designs:
        passed = True
        for k in d['filter_info'].keys():
            #if k == 'fragment_analysis': continue ###DEBUG

            passed = passed and d['filter_info'][k]['pass']

        if passed:
            passed_designs.append(d)

    print 'Number of designs that pass filters = ', len(passed_designs)

    return passed_designs[:max_pass]

def plot_one_filter_score(design_list, filter_name, score_name):
    '''Plot the distribution of a filter score.'''
    data = [d['filter_info'][filter_name][score_name] for d in design_list]
    return PPSD.plot.plot_histogram(data, '_'.join([filter_name, score_name]), show_plot=False) 

def plot_filter_scores(input_path, save_figures=False):
    '''Plot distributions of all scores.'''
    figure_path = os.path.join(input_path, 'figures')
    if not os.path.exists(figure_path):
        os.mkdir(figure_path)
    
    design_list = load_designs(input_path)

    keys = [('SSPrediction', 'score', 0.75, 1),
            ('PackStat', 'score', 0.6, 1),
            ('Holes', 'score', 2, 0),
            ('SSShapeComplementarity', 'score', 0.6, 1),
            ('BuriedUnsatHbonds', 'score', 20, 0),
            ('fragment_analysis', 'worst_crmsd', 1, 0)]
   
    for key in keys:
        plt = plot_one_filter_score(design_list, key[0], key[1])
        PPSD.plot.plot_rectangular_box(key[2], plt.xlim()[key[3]], plt.ylim()[0], plt.ylim()[1],
                savefig_path=os.path.join(figure_path, key[0] + '_distribution.png') if save_figures else None)

def plot_fragment_quality_each_position(input_path, savefig=False):
    '''Plot the fragment quality at each position.'''
    crmsds_list = []

    for d in os.listdir(input_path):
        fragment_discribing_file = os.path.join(input_path, d, 'frags.fsc.200.9mers')
        
        if os.path.exists(fragment_discribing_file):
            crmsds_list.append(
                    PPSD.fragment_quality_analysis.FragmentQualityAnalyzer.get_position_crmsd(
                    fragment_discribing_file))
    
    X = []
    Y = []
    for crmsds in crmsds_list:
        X += list(range(1, len(crmsds) + 1))
        Y += crmsds

    import matplotlib.pyplot as plt

    plt.scatter(X, Y)
    plt.ylim(0, plt.ylim()[1])
    plt.axhline(1, color='r')
    
    plt.xlabel('Sequence position')
    plt.ylabel('CRMSD')
   
    if savefig:
        if not os.path.exists(os.path.join(data_path, 'figures')):
            os.mkdir(os.path.join(data_path, 'figures'))
        
        plt.savefig(os.path.join(data_path, 'figures', 'fragment_quality_each_position.png'))
    else:
        plt.show()


if __name__ == '__main__':
    
    data_path = sys.argv[1]
    
    num_jobs = 1
    job_id = 0
    
    if len(sys.argv) > 3:
        num_jobs = int(sys.argv[2])
        job_id = int(sys.argv[3]) - 1

    pyrosetta.init()

    ###DEBUG

    #pdb_file = 'data/test_assemble/0/assembled.pdb' 
    #pdb_file = '/home/xingjie/DataBases/PPIAffinityDatabase/Affinity_Benchmark/1EMV_r_u.pdb' 
    #pose = rosetta.core.pose.Pose()
    #rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    #filter_dict = {}

    ##ss_prediction_filter(pose, filter_dict)
    ##pack_stat_filter(pose, filter_dict)
    ##holes_filter(pose, filter_dict)
    ##helix_shape_complementarity_filter(pose, filter_dict)
    ##buried_unsatisfied_hbond_filter(pose, filter_dict)

    #print filter_dict

    ####DEBUG

    #filter_designs(data_path, num_jobs, job_id)
    
    plot_filter_scores(data_path)

    #print [(d['id'], d['task_info']['score']) for d in select_designs(data_path, 1000)]

    #plot_fragment_quality_each_position(data_path, savefig=True)
