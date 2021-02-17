#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
import argparse
from peak_stats.reader.read_peaks import Reader
from peak_stats.reader.peaks import Image
from peak_stats.statistics.single_img_stats import PeakPositions, GroupPeakStats


def main():
    parser = argparse.ArgumentParser(
        description="Save statistics for ipalm image peaks file exported from PeakSelector")
    parser.add_argument("-d", "--datafile", help="Ascii file from PeakSelector", type=str, required=True)
    parser.add_argument("-o", "--outdir", help="Where save your plots.")
    parser.add_argument("-c", "--convex", default=False, help="If True plot convex hull in addition to 3d peaks")
    parser.add_argument("-p", "--paint", default=False, help="If true plot 3d peaks in different color for each spot. "
                                                             "Not available for convex hull. Default=False.", type=bool)
    parser.add_argument("-s", "--sigma", default=None, type=int,
                        help="Upper sigma (uncertainty) threshold for X and Y positions of the peak.")
    parser.add_argument("-g", "--group_peaks", type=bool, default=False, help="If True plot group peaks only instead of all peaks.")
    args = parser.parse_args()
    reader = Reader(args.datafile)
    image = Image(reader)
    peaks = PeakPositions(image=image, sigma_threshold=args.sigma, minimize=True)

    print(f"Plotting {args.datafile}")
    if args.group_peaks:
        group_peaks = GroupPeakStats(image=image)
        if args.convex:
            if args.outdir is None:
                outpath = args.datafile[:-4] + "_group_peaks_convex_hull.png"
            else:
                outpath = args.outdir
            group_peaks.plot_group_peaks_convex_hull(outpath=outpath)
        else:
            if args.outdir is None:
                outpath = args.datafile[:-4] + "_group_peaks_3d_plot.png"
            else:
                outpath = args.outdir
            group_peaks.plot_3d_group_peaks(outpath=outpath)
    else:
        if args.convex:
            if args.outdir is None:
                if args.sigma:
                    outpath = args.datafile[:-4] + f"sigma{args.sigma}_convex_hull.png"
                else:
                    outpath = args.datafile[:-4] + "_convex_hull.png"
            else:
                outpath = args.outdir
            peaks.plot_convex_hull(title="Convex hull", outpath=outpath)
        else:
            if args.outdir is None:
                if args.sigma:
                    outpath = args.datafile[:-4] + f"sigma{args.sigma}_3d_peaks.png"
                else:
                    outpath = args.datafile[:-4] + "_3d_peaks.png"
            else:
                outpath = args.outdir
            peaks.plot_peak_positions(title="Peaks 3D plot", outpath=outpath)



if __name__ == '__main__':
    main()
