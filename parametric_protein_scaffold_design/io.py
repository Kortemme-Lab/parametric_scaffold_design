import os
import io
from datetime import timedelta

from flufl.lock import Lock


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
