import numpy as np

def RMSD(points1, poinsts2):
    '''Calcualte RMSD between two lists of numpy points.'''
    diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
    return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))
