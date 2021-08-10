#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
import argparse

from peak_stats.reader.peaks import Image
from peak_stats.reader.read_peaks import Reader
from peak_stats.export_pdb.ascii_2_pdb import Points
from peak_stats.statistics.single_img_stats import GroupPeakStats
from os import path


def main():
    parser = argparse.ArgumentParser(
        description="Save statistics for ipalm image peaks file exported from PeakSelector")
    parser.add_argument("-d", "--datafile", help="Ascii file from PeakSelector", type=str, required=True)
    parser.add_argument("-o", "--output", help="File to save results", default=None, type=str)
    parser.add_argument("-m", "--mode",
                        help="Mode. Which peaks to save.  "
                             " - 'peaks' (saves all peaks to pdb file), 'group' (save group peaks), default='group'",
                        default='group')
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
        img_points.save_pdb(filename=outname + "_all.pdb", points=img_points.points)
    elif args.mode == "group":
        group_peaks = GroupPeakStats(image=image)
        img_points.save_pdb(filename=outname + "_group_peaks.pdb", points=group_peaks.positions)
        pass

    else:
        raise Exception("Wrong mode. Choose from: group or peaks.")


if __name__ == '__main__':
    main()
