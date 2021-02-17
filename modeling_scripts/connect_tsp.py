#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""

# todo - delete spanning tree and emish modeling_operations
import argparse
from modeling_scripts.modeling_operations.create_structure import load_points, save_pdb
from modeling_scripts.modeling_operations.create_structure import run_image_tsp
from scipy.spatial.distance import squareform, pdist


def main():
    parser = argparse.ArgumentParser(
        description="From unorganized set of points (--pdb) identified from density 3D image (--image) create a polymer"
                    " like structure fitting given 3D image using greedy TSP solution.")
    parser.add_argument('-p', '--pdb', help='PDB file with unorganized points', type=str, required=True)
    parser.add_argument('-o', '--optimization_steps', type=int, default=10,
                        help="Maximum number of optimization run on TSP path. default=10")
    args = parser.parse_args()

    if args.pdb.endswith(".pdb"):
        in_points = load_points(infile=args.pdb)
        graph = squareform(pdist(in_points))
        if args.method == "greedy":
            arranged_points = run_image_tsp(graph=graph, points=in_points, endpoints=None,
                                            optimization_steps=args.optimization_steps)
            save_pdb(points=arranged_points, filename=args.pdb[:-4] + "_tsp.pdb", connect=True)

    else:
        raise Exception("Points are required to be in .pdb format.")


if __name__ == '__main__':
    main()
