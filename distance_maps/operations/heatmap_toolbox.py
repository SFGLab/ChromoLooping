#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Set of functions for drawing heatmaps.
Autor: Michał Kadlof <m.kadlof@gmail.com>
"""

import matplotlib.pyplot as plt
from matplotlib import colors
from modeling_scripts.modeling_operations.create_structure import load_points, save_pdb
import numpy as np
import math


# todo change description
# todo improve numpy implementation
def distance_matrix(pdb_file, save=None):
    '''
    read points from pdd file and creat distance matrix between those points
    :param pdb_file: pdb file with points
    :return: distance matrix
    '''
    if not pdb_file.endswith('.pdb'):
        return 'Wrong file'
    points = load_points(pdb_file)
    n = len(points)
    distances = np.zeros((n, n))
    for i in range(len(points)):
        for j in range(len(points)):
            x1, y1, z1 = points[i][0], points[i][1], points[i][2]
            x2, y2, z2 = points[j][0], points[j][1], points[j][2]
            distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
            distances[i][j] = distance
            distances[j][i] = distance
    if save:
        np.savetxt(save, distances, fmt='%d')
        print("File {} saved...".format(save))
    else:
        return distances


# todo work on options and clean this + description
def show(heatmap, gamma=1, draw_line1=None, drawline2=None, title="No title", interpolation="nearest",
         out_file_name=None, cmap='afmhot', vmax=None):
    """Plots beautiful heatmap.

       Args:
           heatmap (np.array): 2d numpy array
           gamma (float): scalling factor in power law used in color mapping. Default: 1 (linear)
           genomic_postion_start (int): adjust axis scales
           genomic_postion_end (int): adjust axis scales

    """
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    if vmax:
        cax = ax1.imshow(heatmap, interpolation="nearest", cmap=cmap, norm=colors.PowerNorm(gamma=gamma), aspect=1.0,
                         vmax=vmax)
    else:
        cax = ax1.imshow(heatmap, interpolation="nearest", cmap=cmap, norm=colors.PowerNorm(gamma=gamma), aspect=1.0)
    fig.colorbar(cax)
    plt.title(title)
    if draw_line1:
        plt.axvline(x=1000, color="yellow")
    if drawline2:
        plt.axhline(y=2331, color="cyan")
    ax1.set_xlabel('model bead index')
    ax1.set_ylabel('model bead index')
    # chromosome positions i ccd scale (whole genome)
    # ticks = [189,368, 534,653,786,923,1053,1164,1258,1368,1476,1592,1662,1733,1805,1872,1955,2009,2066,2112,2134,2166] #just for whole genome!!!!!
    # plt.setp(ax1.get_xticklabels(), rotation='vertical', fontsize=6)
    # plt.setp(ax1.get_yticklabels(), fontsize=6)
    # plt.xticks(ticks)
    # plt.yticks(ticks)
    # plt.grid(color='w', linestyle='-',linewidth=.25)

    fig.subplots_adjust(bottom=0.2, left=0.2)
    '''
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    ax2.spines["bottom"].set_position(("axes", -0.1))
    ax2.set_xticklabels(tick_values)
    ax2.set_xticks(tick_locations)
    ax2.set_xlabel('genomic postion [Mbp]')

    ax3.yaxis.set_ticks_position("left")
    ax3.yaxis.set_label_position("left")
    ax3.spines["left"].set_position(("axes", -0.1))
    ax3.set_yticklabels(tick_values)
    ax3.set_yticks(tick_locations)
    ax3.set_ylabel('genomic postion [Mbp]')
    '''
    if out_file_name:
        # fig.set_size_inches(18.78, 15, forward=True)
        plt.savefig(out_file_name, bbox_inches='tight', dpi=300)
    else:
        plt.show()
    plt.close()
    # return plt
