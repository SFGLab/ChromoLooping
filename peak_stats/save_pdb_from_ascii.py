#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
import argparse

from peak_stats.reader.peaks import Image
from peak_stats.reader.read_peaks import Reader
from peak_stats.export_pdb.ascii_2_pdb import Points
from os import path


def main():
    parser = argparse.ArgumentParser(
        description="Save statistics for ipalm image peaks file exported from PeakSelector")
    parser.add_argument("-d", "--datafile", help="Ascii file from PeakSelector", type=str, required=True)
    parser.add_argument("-o", "--output", help="File to save results", default=None, type=str)
    parser.add_argument("-m", "--mode",
                        help="Mode. Which peaks to save. All or a representant od a group. peaks (saves all peaks to "
                             "pdb file), groups(save peaks of each group in separate file), 'center' (save cetroids of "
                             "all groups) or 'all' (save all of the options)", default="peaks")
    parser.add_argument("-s", "--sigma", default=None,
                        help="Upper sigma (uncertainty) threshold for X and Y positions of the peak.", type=float)
    args = parser.parse_args()
    reader = Reader(args.datafile)
    image = Image(reader)
    img_points = Points(image=image, sigma_threshold=args.sigma)
    if not args.output:
        outname = path.join(path.dirname(args.datafile), path.basename(args.datafile).split(".")[0])
    else:
        outname = args.output
    if args.mode == "peaks":
        # all peaks
        img_points.save_pdb(filename=outname + "_all.pdb", points=img_points.points)
    elif args.mode == "center":
        # centroinds for each group
        centroids = img_points.parse_image()
        img_points.save_pdb(filename=outname + "_centroids.pdb", points=centroids)
    elif args.mode == "groups":
        # each group in different pdb file
        img_points.save_spots(outname=outname, sigma_threshold=args.sigma)
    elif args.mode == "all":
        img_points.save_pdb(filename=outname + "_all.pdb", points=img_points.points)
        centroids = img_points.parse_image(sigma_threshold=args.sigma)
        img_points.save_pdb(filename=outname + "_centroids.pdb", points=centroids)
        img_points.save_spots(outname=outname, sigma_threshold=args.sigma)
    elif args.mode == "oligos":
        # todo dodac opcje zapisywania group peaks
        pass

    else:
        raise Exception("Wrong mode. Choose from: all, center, groups and peaks.")


if __name__ == '__main__':
    main()
