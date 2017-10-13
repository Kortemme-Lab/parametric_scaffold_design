'''Functions to analyze sequences.'''
import subprocess


def simple_sequence_identity(sequence1, sequence2):
    '''Calculate the sequence identity of two sequences with the same length.'''
    num_identity = len([i for i in range(len(sequence1)) if sequence1[i] == sequence2[i]])
    return num_identity * 1.0 / len(sequence1)


def combine_fasta_files(input_files, combined_file):
    '''Combine a list of fasta files into a single file'''
    with open(combined_file, 'w') as cf:
        for input_file in input_files:
            with open(input_file, 'r') as i_f:
                cf.write(i_f.read() + '\n')

def make_sequence_logo(sequence_file, image_file):
    '''Make sequence logos for a given (fasta) sequence file using the weblogo application.'''
    cmd = ['weblogo', '-c', 'chemistry', '--format', 'PNG', '--fin', sequence_file, '--fout', image_file]
    
    subprocess.check_call(cmd) 
