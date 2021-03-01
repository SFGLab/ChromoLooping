#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
functions to perform spline interpolation
"""

from scipy import interpolate
import numpy as np
from .distance import distance3d as d
from .distance import point_between
import random


# todo - removing points is not very elegant. Is there maybe a better way
def spline(points, ref_dist=1.0, max_beads=None):
    """Performs spline interpolation out of list of points."""
    points = np.asarray(points)
    points = points.transpose()
    tck, u = interpolate.splprep(points, s=0)
    new = interpolate.splev(np.linspace(0, 1, 100000), tck)
    new_l = list(zip(new[0], new[1], new[2]))
    p = [new_l[0]]
    n = len(new_l)
    previous_point = new_l[0]
    distance_to_last_bead = 0
    for i in range(n):
        current_point = new_l[i]
        loc_dist = d(previous_point, current_point)
        distance_to_last_bead += loc_dist
        if distance_to_last_bead > ref_dist:
            proportion = (distance_to_last_bead - ref_dist) / loc_dist
            p.append(point_between(previous_point, current_point, prop=proportion))
            distance_to_last_bead = distance_to_last_bead - ref_dist
        previous_point = current_point
    if max_beads:
        while len(p) > max_beads:
            to_remove = random.choice(p)
            if to_remove != p[0] and to_remove != p[-1]:
                p.remove(to_remove)
    return p

def refdist(points, beads):
    """Calculate polymer length"""
    structure_length = 0
    for i in range(len(points) - 1):
        point_1 = points[i]
        point_2 = points[i + 1]
        dist = d(point_1, point_2)
        structure_length += dist
    ref_dist = structure_length / float(beads)
    return ref_dist