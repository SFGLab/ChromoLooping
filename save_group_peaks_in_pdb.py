#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
import argparse

from peak_stats.reader.peaks import Image, GroupPeak
from peak_stats.reader.read_peaks import Reader
from peak_stats.statistics.single_img_stats import GroupPeakStats
from peak_stats.export_pdb.ascii_2_pdb import Points
from os import path


def main():
    parser = argparse.ArgumentParser(
        description="Save statistics for ipalm image peaks file exported from PeakSelector")
    parser.add_argument("-d", "--datafile", help="Ascii file from PeakSelector", type=str, required=True)
    parser.add_argument("-o", "--output", help="File to save results", default=None, type=str)
    parser.add_argument("-s", "--sigma", default=None,
                        help="Upper sigma (uncertainty) threshold for X and Y positions of the peak.", type=float)
    args = parser.parse_args()
    reader = Reader(args.datafile)
    image = Image(reader)
    group_peaks = GroupPeakStats(image=image)
    points = Points(image=image, sigma_threshold=args.sigma)
    if args.output:
        points.save_pdb(args.output, points=group_peaks.positions)
    else:
        output = args.datafile[:-4] + "_group_peaks.pdb"
        points.save_pdb(filename=output, points=group_peaks.positions)


if __name__ == '__main__':
    main()
