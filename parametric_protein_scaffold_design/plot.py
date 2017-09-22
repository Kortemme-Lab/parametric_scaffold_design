import os

import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt


def plot_histogram(data, xlabel, x_min=None, x_max=None, savefig_path=None, savefig_type='png', show_plot=True):
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

    if savefig_path:
        plt.savefig(os.path.join(savefig_path, '_'.join(xlabel.split()) + '.' + savefig_type))
    elif show_plot:
        plt.show()

    return plt

def plot_rectangular_box(left, right, bottom, top, color='r', linewidth=3.0, savefig_path=None, show_plot=True):
    '''Plot a rectangular box'''

    plt.plot([left, left, right, right, left], [bottom, top, top, bottom, bottom], color=color, linewidth=linewidth)

    if savefig_path:
        plt.savefig(savefig_path)
    elif show_plot:
        plt.show()
