#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""

import argparse
from distance_maps.operations.heatmap_toolbox import distance_matrix, show
from os import path, listdir
import numpy as np


def dir_exist(path_to_dir):
    if not path.isdir(path_to_dir):
        raise argparse.ArgumentTypeError("Wrong path")
    return path_to_dir


def main():
    parser = argparse.ArgumentParser(description="opis")
    parser.add_argument('pdb', type=str, help='PDB file or path to multiple pdb files')
    parser.add_argument('-s', '--save', default=None, help='Path to save your output files. If not provided files are '
                                                           'saved in input direction', type=dir_exist)
    parser.add_argument('-i', '--image', type=bool, help='If True - save .png files for heatmaps', default=False)
    parser.add_argument('-c', '--cmap', type=str, default='seismic_r',
                        help='Python color map for heatmaps. Only if you '
                             'want to save png files.')
    parser.add_argument('-r', '--reverse', type=bool, default=False,
                        help="Save two distance maps, one from input model, "
                             "other one from reverse model, useful for image "
                             "driven modelling, when start and end position is unknown")
    parser.add_argument("-t", "--title", default="Distance map", help="Title of the ploted distance map. Used only if "
                                                                      "-i is used. Default is 'Distance map'. ")
    args = parser.parse_args()
    outfiles = []
    if path.isdir(args.pdb):
        files = listdir(args.pdb)
        for i in files:
            if i.endswith('.pdb'):
                infile = path.join(args.pdb, i)
                if args.save:
                    outfile = path.join(args.save, i[:-4] + '.heat')
                else:
                    outfile = path.join(args.pdb, i[:-4] + '.heat')
                outfiles.append(outfile)
                distance_matrix(infile, outfile)
    else:
        infile = args.pdb
        if args.save:
            outfile = path.join(args.save, args.pdb[:-4] + '.heat')
        else:
            outfile = path.join(args.pdb, args.pdb[:-4] + '.heat')
        outfiles.append(outfile)
        distance_matrix(infile, outfile)
    if args.image:
        for i in outfiles:
            heatmap = np.loadtxt(i)
            image_outfile = f"{i[:-5]}_{args.cmap}.png"
            title = 'Distance map'
            show(heatmap, gamma=1, out_file_name=image_outfile, title=title, cmap=args.cmap)


if __name__ == '__main__':
    main()
