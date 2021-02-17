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


def spline(points, ref_dist=1.0, max_beads=None):
    """Performs spline interpolation out of list of points

    Args:
        points (list of three floats tuples) : points that will be used as spline nodes
        ref_dist (float) : desired distances between beads

    Returns:
        points (list of three floats tuples) : points after interpolation
    """
    points = np.asarray(points)
    points = points.transpose()
    tck, u = interpolate.splprep(points, s=0)
    new = interpolate.splev(np.linspace(0, 1, 100000), tck)
    bla = list(zip(new[0], new[1], new[2]))
    p = [bla[0]]
    n = len(bla)
    previous_point = bla[0]
    distance_to_last_bead = 0
    for i in range(n):
        current_point = bla[i]
        loc_dist = d(previous_point, current_point)
        distance_to_last_bead += loc_dist
        if distance_to_last_bead > ref_dist:
            proportion = (distance_to_last_bead - ref_dist)/loc_dist
            p.append(point_between(previous_point, current_point, prop=proportion))
            distance_to_last_bead = distance_to_last_bead - ref_dist
        previous_point = current_point
    if max_beads:
        while len(p) > max_beads:
            to_remove = random.choice(p)
            if to_remove != p[0] and to_remove != p[-1]:
                p.remove(to_remove)
    return p

