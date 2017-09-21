import os
import io
import json
from datetime import timedelta

from flufl.lock import Lock

import pyrosetta
import pyrosetta.rosetta as rosetta


def safe_append(file_name, string):
    '''Append a string to a file. This function
    is thread safe.
    '''
    # Make a lock

    lock_name = os.path.join(file_name + '.lock')
    lock = Lock(lock_name)
    lock.lifetime = timedelta(minutes=10)

    # Write to the result file with a lock
    
    with lock:
        with open(file_name, 'a+') as f:
            f.write(string)

def ll_to_tsv(file_name, ll):
    '''Save a list of list into a tsv file.'''
    with open(file_name, 'w') as f:
        for l in ll:
            f.write('\t'.join([repr(x) for x in l]) + '\n')

def tsv_to_ll(file_name, cast_funs):
    '''Read a tsv file into a list of file. A list of type casting
    functions will be applied to each line of the file.'''
    ll = []

    with open(file_name, 'r') as f:
        for line in f.readlines():
            ll.append([])
            sline = line.split()

            for i in range(len(cast_funs)):
                ll[-1].append(cast_funs[i](sline[i]))

    return ll

def write_rigid_body_transformation_to_file(T, fout):
    '''Write a rigid body transformation T into a json file.'''
    with open(fout, 'w') as f:
        json.dump([[T.xx(), T.yx(), T.zx(), T.x()],
            [T.xy(), T.yy(), T.zy(), T.y()],
            [T.xz(), T.yz(), T.zz(), T.z()]], f)
    
def load_rigid_body_transformation_from_file(fin):
    '''Load a rigid body transformation from a json file.'''
    T = rosetta.numeric.xyzTransform_double_t()
  
    with open(fin, 'r') as f:
        l = json.load(f)
        
        T.R.xx = l[0][0]
        T.R.yx = l[0][1]
        T.R.zx = l[0][2]
        T.t.x = l[0][3]

        T.R.xy = l[1][0]
        T.R.yy = l[1][1]
        T.R.zy = l[1][2]
        T.t.y = l[1][3]

        T.R.xz = l[2][0]
        T.R.yz = l[2][1]
        T.R.zz = l[2][2]
        T.t.z = l[2][3]

    return T

def sequence_to_fasta_file(fasta_path, title, sequence):
    '''Save a sequence into a fasta file.'''
    with open(fasta_path, 'w') as f:
        f.write('>' + title + '\n')
        
        start = 0
        end = min(80, len(sequence))
        f.write(sequence[start:end] + '\n')

        while end < len(sequence):
            start += 80
            end = min(start + 80, len(sequence))
            f.write(sequence[start:end] + '\n')

