'''Functions to analyze sequences.'''
import subprocess
import matplotlib.pyplot as plt



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

def make_sequence_logo(sequence_file, image_file, output_format='png_print'):
    '''Make sequence logos for a given (fasta) sequence file using the weblogo application.'''
    cmd = ['weblogo', '-c', 'chemistry', '--format', output_format, '--fin', sequence_file, '--fout', image_file]
    
    subprocess.check_call(cmd)

def filter_sequences_by_identity(sequences, threshold):
    '''Get a subset of given sequences in which no sequence identity is above a threshold.
    Sequences are readin one by one and if the identity between a new sequence and any sequence
    that is already saved in the new set is above the threshold, the new sequence will not be
    saved.
    '''
    non_redundant_sequences = []
    for seq in sequences:
        save = True
        for saved_seq in non_redundant_sequences:
            if simple_sequence_identity(seq, saved_seq) > threshold:
                save = False
                break

        if save:
            non_redundant_sequences.append(seq)

    return non_redundant_sequences

def plot_sequence_identities(sequences):
    '''Make a plot of the number of sequences under different identity cutoffs.
    Note that the order of the sequences will affect the result.
    '''
    cut_offs = [0.1 * i for i in range(11)]
    num_sequences = []

    for cut_off in cut_offs:
        num_sequences.append(len(filter_sequences_by_identity(sequences, cut_off)))

    plt.plot(cut_offs, num_sequences)
    plt.xlim(0, 1)
    plt.ylim(0, len(sequences))
    plt.xlabel('sequence identity')
    plt.ylabel('number of sequences')

    plt.show()
