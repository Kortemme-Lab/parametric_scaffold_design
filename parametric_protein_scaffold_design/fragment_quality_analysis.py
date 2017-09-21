import os
import subprocess


class FragmentQualityAnalyzer:
    '''Analyze fragment quality of a strcuture.'''
    def __init__(self, runpsipred_single_path, fragment_picker_cmd, vall_path, weights_file):
        self.runpsipred_single_path_abs = os.path.abspath(runpsipred_single_path)
        self.fragment_picker_cmd = fragment_picker_cmd
        self.vall_path_abs = os.path.abspath(vall_path)
        self.weights_file_abs = os.path.abspath(weights_file)
    
    def pick_fragments(self, input_pdb, input_fasta, working_dir):
        '''Pick fragments of a structure.
        Return the name of the fragment describing file.
        '''
        cwd = os.getcwd()
        input_pdb_abs = os.path.abspath(input_pdb)
        input_fasta_abs = os.path.abspath(input_fasta)
    
        os.chdir(working_dir)
    
        # Use PSIPRED to generate secondary structure predictions
    
        subprocess.check_call([self.runpsipred_single_path_abs, input_fasta_abs])
        ss2_file = os.path.basename(input_fasta)[:-5] + 'ss2'
    
        # Do fragment picking
   
        pick_cmd = [self.fragment_picker_cmd,
                    '-in::file::vall', self.vall_path_abs,
                    '-in::file::fasta', input_fasta_abs,
                    '-in::file::s', input_pdb_abs,
                    '-frags::ss_pred', ss2_file, 'predA',
                    '-frags::scoring::config', self.weights_file_abs,
                    '-frags::frag_sizes', '9',
                    '-frags::n_candidates', '200',
                    '-frags::n_frags', '200',
                    '-out::file::frag_prefix', 'frags',
                    '-frags::describe_fragments', 'frags.fsc']

        subprocess.check_call(pick_cmd)


        # Go back to the previous working directory
    
        os.chdir(cwd)

        return os.path.join(working_dir, 'frags.fsc.200.9mers')

    def get_position_crmsd(self, fragment_discribing_file):
        '''Return a list of best fragment crmsd at each position.'''
        best_crmsds = []
        
        with open(fragment_discribing_file, 'r') as fdf:
            for line in fdf.readlines():
                
                if line.startswith('#'):
                    best_crmsds.append(float('inf'))
                    continue
                
                crmsd = float(line.split()[-3])

                if crmsd < best_crmsds[-1]:
                    best_crmsds[-1] = crmsd
   
        return best_crmsds

if __name__ == '__main__':

    fqa = FragmentQualityAnalyzer('./dependencies/dependencies/psipred/runpsipred_single', 'fragment_picker.linuxclangrelease', 'database/fragment_quality_analysis/small.vall.gz', 'database/fragment_quality_analysis/simple.wghts')

    #fdf = fqa.pick_fragments('data/fragment_quality_analysis/assembled.pdb', 'data/fragment_quality_analysis/asseA.fasta', 'data/fragment_quality_analysis')

    fdf = 'data/fragment_quality_analysis/frags.fsc.200.9mers'
    print fqa.get_position_crmsd(fdf)

