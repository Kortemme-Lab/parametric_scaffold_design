import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt


def plot_histogram(data, xlabel, x_min=None, x_max=None):
    '''Plot a histogram.'''
    if x_min is None:
        x_min = min(data)
    
    if x_max is None:
        x_max = max(data)
    
    step = (x_max - x_min) / 100.0
    hist, bin_edges = np.histogram(data, bins=np.arange(x_min, x_max, step))

    plt.clf()
    plt.bar(bin_edges[0:-1], hist, width=step, edgecolor='black')
    plt.ylabel('Number of data')
    plt.xlabel(xlabel)

    plt.show()
