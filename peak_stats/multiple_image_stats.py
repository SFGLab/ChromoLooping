#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""

import argparse
import glob
from peak_stats.statistics.multiple_img_stats import ImageDict


def main():
    parser = argparse.ArgumentParser(
        description="Calculate statistics for multiple images (ASCII files) from PeakSelector.")
    parser.add_argument("-p", "--path",
                        help="Directory to search using glob. Example: /mnt/raid/images/*/*/*[0-9].txt",
                        type=str, required=True)
    parser.add_argument("-o", "--output", help="Directory to save statistics in.", default=None, type=str)
    parser.add_argument("-s", "--save_plots", help="If True save histograms in given directory", default=False,
                        type=bool)
    args = parser.parse_args()

    file_list = glob.glob(args.path)
    im_dict = ImageDict.create_from_files(inpath=file_list)
    im_dict.save_stats(output=args.output)
    if args.save_plots:
        im_dict.save_plots(output=args.output)


if __name__ == '__main__':
    main()
