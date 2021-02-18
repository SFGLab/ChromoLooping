#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""

from peak_stats.reader.read_peaks import Reader
from peak_stats.reader.peaks import Image
from peak_stats.statistics.single_img_stats import ImgStats
import argparse


# todo - write in description which statistics are saved in text files
def main():
    parser = argparse.ArgumentParser(
        description="Save statistics (in text file) for ipalm image peaks file exported from PeakSelector.")
    parser.add_argument("-d", "--datafile", help="Ascii file from PeakSelector", type=str, required=True)
    parser.add_argument("-o", "--output", help="File to save results", default=None, type=str)
    parser.add_argument("-s", "--sigma", default=None,
                        help="Upper sigma (uncertainty) threshold for X and Y positions of the peak.")
    args = parser.parse_args()
    reader = Reader(args.datafile)
    image = Image(reader)
    stats = ImgStats(image=image)
    if args.output:
        stats.save_statistics(output=args.output, image=image, sigma_threshold=args.sigma)
    else:
        output = args.datafile[:-4] + "_stats.txt"
        stats.save_statistics(output=output, image=image, sigma_threshold=args.sigma)


if __name__ == '__main__':
    main()
