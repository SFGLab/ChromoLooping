#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 15:41:22 2017

@author: zofia
"""
import argparse
from modeling_scripts.modeling_operations.create_structure import load_points, save_pdb
from modeling_scripts.smoothing_operations.spline import spline
from modeling_scripts.smoothing_operations.distance import distance3d as d


def main():
    parser = argparse.ArgumentParser(description="Run spline interpolation on polymer model in .pdb format.")
    parser.add_argument("-i", '--infile', type=str, help='model to smooth')
    parser.add_argument('-b', '--beads', type=int, help='number of beads in final model')
    args = parser.parse_args()
    out = args.infile[:-4] + f'_smooth{args.beads}.pdb'
    points = load_points(args.infile)
    ref_dist = refdist(points, args.beads)
    new_points = spline(points, ref_dist, max_beads=args.beads)
    save_pdb(new_points, out)


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


if __name__ == '__main__':
    main()
