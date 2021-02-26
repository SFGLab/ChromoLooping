#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
from tsp_solver.greedy import solve_tsp
import numpy as np
import re


def run_image_tsp(graph, points, endpoints=None, optimization_steps=10):
    """Run greedy TSP on input graph, and with optionally fixed start and endpoint ex. (1,3)"""
    if endpoints:
        path = solve_tsp(graph, optim_steps=optimization_steps, endpoints=endpoints)
    else:
        path = solve_tsp(graph)
    arranged_points = [0] * len(path)
    for i in range(len(path)):
        arranged_points[i] = points[path[i]]
    return arranged_points


def load_points(infile):
    """Load points form .PDB"""
    atoms = [i.strip() for i in open(infile) if re.search('^(ATOM|HETATM)', i)]
    points = []
    for i in atoms:
        x = float(i[30:38])
        y = float(i[38:46])
        z = float(i[46:54])
        points.append((x, y, z))
    return np.array(points)


def save_pdb(points, filename, verbose=True, connect=True):
    """Save points in .PDB format."""
    atoms = ''
    n = len(points)
    for i in range(n):
        x = points[i][0]
        y = points[i][1]
        try:
            z = points[i][2]
        except IndexError:
            z = 0.0
        atoms += (
            '{0:6}{1:>5}  {2:3}{3:}{4:3} {5:}{6:>4}{7:}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}{13:>12}\n'.format(
                "ATOM", i + 1, 'B', ' ', 'BEA', 'A', i + 1, ' ', x, y, z, 0, 0, 'C'))

    connects = ''
    if connect:
        if n != 1:
            connects = 'CONECT    1    2\n'
            for i in range(2, n):
                connects += 'CONECT{:>5}{:>5}{:>5}\n'.format(i, i - 1, i + 1)
            connects += 'CONECT{:>5}{:>5}\n'.format(n, n - 1)

    pdb_file = atoms + connects
    open(filename, 'w').write(pdb_file)
    if verbose:
        print("File {} saved...".format(filename))
    return filename
