import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt


def plot_histogram(data, xlabel):
    '''Plot a histogram.'''
    step = (max(data) - min(data)) / 100
    hist, bin_edges = np.histogram(data, bins=np.arange(min(data), max(data), step))

    plt.clf()
    plt.bar(bin_edges[0:-1], hist, width=step, edgecolor='black')
    plt.ylabel('Number of data')
    plt.xlabel(xlabel)

    plt.show()
