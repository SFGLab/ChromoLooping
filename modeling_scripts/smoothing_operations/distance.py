#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module with some usefull functions for manipulating points"""

from math import sqrt


def distance3d(point_1, point_2):
    """Calculates euclidean distance in R^3 between two points

    Args:
        point_1 (three floats tuple) : point 1 coordinates
        point_2 (three floats tuple) : point 2 coordinates

    Returns:
        float: euclidean distance between point 1 and point 2
    """
    return sqrt((point_1[0]-point_2[0])**2+(point_1[1]-point_2[1])**2+(point_1[2]-point_2[2])**2)


def distance2d(point_1, point_2):
    """Calculates euclidean distance in R^2 between two points

    Args:
        point_1 (two floats tuple) : point 1 coordinates
        point_2 (two floats tuple) : point 2 coordinates

    Returns:
        float: euclidean distance between point 1 and point 2
    """
    return sqrt((point_1[0]-point_2[0])**2+(point_1[1]-point_2[1])**2)


def point_between(point_1, point_2, prop=0.5):
    """Calculates coordiantes of in the middle of two points in R^3.

    Args:
        point_1 (three floats tuple) : point 1 coordinates
        point_2 (three floats tuple) : point 2 coordinates

    Returns:
        (three floats tuple) : point in the middle
    """
    x_1, y_1, z_1 = point_1
    x_2, y_2, z_2 = point_2
    return ( x_1*(1-prop) + x_2 * prop, y_1 * (1-prop) + y_2 * prop, z_1 * (1-prop) + z_2 * prop )

