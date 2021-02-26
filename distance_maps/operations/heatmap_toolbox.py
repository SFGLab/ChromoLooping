#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""

import matplotlib.pyplot as plt
from matplotlib import colors
from modeling_scripts.modeling_operations.create_structure import load_points, save_pdb
import numpy as np
import math


def distance_matrix(pdb_file, save=None):
    """Read points from pdb file and create distance matrix between them"""
    if not pdb_file.endswith('.pdb'):
        raise Exception('Wrong file format.')
    points = load_points(pdb_file)
    n = len(points)
    distances = np.zeros((n, n))
    for i in range(len(points)):
        for j in range(len(points)):
            x1, y1, z1 = points[i][0], points[i][1], points[i][2]
            x2, y2, z2 = points[j][0], points[j][1], points[j][2]
            distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
            distances[i, j] = distance
            distances[j, i] = distance
    if save:
        np.savetxt(save, distances, fmt='%d')
        print("File {} saved...".format(save))
    else:
        return distances


def show(heatmap, gamma=1, draw_x_line=None, draw_y_line=None, title=None, out_file_name=None, cmap='afmhot', ylabel=None, xlabel=None):
    """Plots beautiful heatmap."""
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cax = ax1.imshow(heatmap, interpolation="nearest", cmap=cmap, norm=colors.PowerNorm(gamma=gamma), aspect=1.0)
    fig.colorbar(cax)
    if title:
        plt.title(title)
    if draw_x_line:
        plt.axvline(x=draw_x_line, color="yellow")
    if draw_y_line:
        plt.axhline(y=draw_y_line, color="cyan")
    if ylabel:
        ax1.set_ylabel(ylabel)
    else:
        ax1.set_ylabel("beads")
    if xlabel:
        ax1.set_xlabel(xlabel)
    else:
        ax1.set_xlabel("beads")
    fig.subplots_adjust(bottom=0.2, left=0.2)
    if out_file_name:
        plt.savefig(out_file_name, bbox_inches='tight', dpi=300)
    else:
        plt.show()
    plt.close()
