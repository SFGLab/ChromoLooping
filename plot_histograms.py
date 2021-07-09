#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
import argparse
from peak_stats.reader.read_peaks import Reader
from peak_stats.reader.peaks import Image
from peak_stats.statistics.single_img_stats import ImgStats
from peak_stats.statistics.single_img_stats import PeakPositions


def main():
    parser = argparse.ArgumentParser(
        description="Plot histograms for sigma parameters for X, Y, Z or average sigma for spots")
    parser.add_argument("-d", "--datafile", help="Ascii file from PeakSelector", type=str, required=True)

    parser.add_argument("-p", "--parameter",
                        help="Which parameter would you like to plot a histogram from."
                             "'x' - distribution of uncertainty of x positions for all peaks"
                             "'y' - distribution of uncertainty of y positions for all peaks"
                             "'z' - distribution of uncertainty of z positions for all peaks"
                             "'peaks' - distribution of number of peaks in groups"
                             "'average_peaks' - distribution of average sigma (X, Y, Z) for all peaks"
                             "'avg' - distribution of average sigma (X, Y, Z) for all groups"
                             "'photons' - distribution of number of photons for all peaks"
                             "'avg_photons' - distribution of average number of photons for all groups"
                             "'all' - saves all above plots in given location"
                             " Default is x.",
                        default='x', type=str)
    parser.add_argument("-o", "--outdir", type=str, help="Outdir to save plots.", default=None)
    args = parser.parse_args()
    reader = Reader(args.datafile)
    image = Image(reader)
    stats = ImgStats(image=image)
    stats.add_groups_sigma(image=image)
    sigmas = ImgStats.add_peak_sigma(image)
    if args.parameter == "x":
        stats.plot_sigma_x_hist()
    elif args.parameter == "y":
        stats.plot_sigma_y_hist()
    elif args.parameter == "z":
        stats.plot_sigma_z_hist()
    elif args.parameter == "avg":
        stats.plot_average_sigma()
    elif args.parameter == "peaks":
        stats.plot_peaks_per_group_histogram(image)
    elif args.parameter == "avg_photons":
        stats.plot_photons_per_group()
    elif args.parameter == "photons":
        stats.plot_photons_per_peaks()
    elif args.parameter == "all":
        stats.plot_sigma_x_hist(save=True, outdir=args.outdir)
        stats.plot_sigma_y_hist(save=True, outdir=args.outdir)
        stats.plot_sigma_z_hist(save=True, outdir=args.outdir)
        stats.plot_average_sigma(save=True, outdir=args.outdir)
        stats.plot_peaks_per_group_histogram(image, save=True, outdir=args.outdir)
        stats.plot_photons_per_group(save=True, outdir=args.outdir)
        stats.plot_photons_per_peaks(save=True, outdir=args.outdir)
    elif args.parameter == "average_peaks":
        if args.outdir is None:
            outdir = args.datafile[:-4]
        else:
            outdir = args.outdir
        stats.plot_average_peak_sigma(sigma_z=sigmas[2], sigma_x=sigmas[0], sigma_y=sigmas[1], save=True, outdir=outdir)
    elif args.parameter == "peaks_sigma":
        if args.outdir is None:
            outdir = args.datafile[:-4]
        else:
            outdir = args.outdir
        stats.plot_peak_sigma_x(x_sigma=sigmas[0], save=True, outdir=outdir)
        stats.plot_peak_sigma_y_(y_sigma=sigmas[1], save=True, outdir=outdir)
        stats.plot_peak_sigma_z(sigma_z=sigmas[2], save=True, outdir=outdir)

    else:
        raise ValueError("Wrong parameter. Choose from 'x', 'y', 'z', 'peaks' or 'avg'")


if __name__ == '__main__':
    main()
