#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""

from collections.abc import MutableMapping
from typing import Iterator, Tuple, Dict
from copy import copy
from os import path

from peak_stats.statistics.single_img_stats import PeakPositions, ImgStats
from peak_stats.reader.read_peaks import Reader
from peak_stats.reader.peaks import Image
import glob
import matplotlib.pyplot as plt

# Histograms global options
FIGSIZE = (30, 20)
LABELFONTSIZE = 25
TITLEFONTSIZE = 30
TICKFONTSIZE = 20
FACECOLOR = "g"
ALPHA = 0.75
DPI = 100

ImageDictValue = Tuple[PeakPositions, ImgStats]


class ImageDict(MutableMapping):

    def __delitem__(self, k: Image) -> None:
        del self.images[k]

    def __getitem__(self, k: Image) -> ImageDictValue:
        return self.images[k]

    def __len__(self) -> int:
        return len(self.images)

    def __iter__(self) -> Iterator[Image]:
        return iter(self.images)

    def __setitem__(self, k: Image, v: ImageDictValue) -> None:
        self.images[k] = v

    def __init__(self, images: Dict[Image, ImageDictValue]):
        self.images = copy(images)

    def plot_histogram(self, data, title: str, xlabel: str, ylabel: str, outname=None):
        """Base to plottng multiple histograms"""
        plt.figure(figsize=FIGSIZE)
        plt.hist(data, 40, facecolor=FACECOLOR, alpha=ALPHA)
        plt.ylabel(ylabel, fontsize=LABELFONTSIZE)
        plt.xlabel(xlabel, fontsize=LABELFONTSIZE)
        plt.title(title, fontsize=TITLEFONTSIZE)
        plt.tick_params(axis='both', which='major', labelsize=TICKFONTSIZE)
        plt.grid(True)
        if outname:
            plt.savefig(outname, dpi=DPI, formay="png")
        else:
            plt.show()
        plt.close()

    def save_stats(self, output):
        """Save some statistics for all images."""
        outfile = path.join(output, "all_images_stats.txt")
        with open(outfile, 'w') as out:
            out.write("Number of images: {}".format(len(self.images)) + "\n")
            out.write("Average number of photons per spot: {}".format(self.avg_photons_per_spot()) + "\n")
            out.write("Average number of peaks per image: {}".format(self.avg_peaks_per_image()) + "\n")
            out.write("Average number of peaks per group: {}".format(self.avg_peaks_per_spot()) + "\n")
            out.write("Average number of groups per image: {}".format(self.avg_groups_per_image()))

    def print_stats(self):
        print("Number of images: {}".format(len(self.images)) + "\n")
        print("Average number of photons per spot: {}".format(self.avg_photons_per_spot())  )
        print("Average number of peaks per image: {}".format(self.avg_peaks_per_image()))
        print("Average number of peaks per group: {}".format(self.avg_peaks_per_spot()))
        print("Average number of groups per image: {}".format(self.avg_groups_per_image()))

    def save_plots(self, output):
        # avg photons per spot
        name1 = path.join(output, "all_images_photons_per_group.png")
        data1 = self.photons_per_spot()
        self.plot_histogram(data=data1, title="Average number of photons per group", xlabel="Number of photons",
                            ylabel="Number of groups", outname=name1)
        # avg peak count
        name2 = path.join(output, "all_images_peaks_per_image.png")
        data2 = self.peaks_per_image()
        self.plot_histogram(data=data2, title="Number of peaks per image", xlabel="Number of peaks",
                            ylabel="Number of images", outname=name2)
        # peaks per spot
        name3 = path.join(output, "all_images_peaks_per_spot.png")
        data3 = self.peaks_per_spot()
        self.plot_histogram(data=data3, title="Peaks per spot", xlabel="Number of peaks", ylabel="Number of spots",
                            outname=name3)

    def avg_photons_per_spot(self):
        photons_per_spot = 0
        for i in self.images.keys():
            photons_per_spot += sum(self.images[i][1].photons_per_spot) / len(self.images[i][1].photons_per_spot)
        return photons_per_spot / len(self)

    def photons_per_spot(self):
        photons_per_spot = []
        for i in self.images.keys():
            photons_per_spot += self.images[i][1].photons_per_spot
        return photons_per_spot

    def peaks_per_image(self):
        peaks_per_image = []
        for i in self.images.keys():
            peaks_per_image.append(self.images[i][1].peak_count)
        return peaks_per_image

    def avg_peaks_per_image(self):
        peaks_per_image = self.peaks_per_image()
        avg = sum(peaks_per_image) / len(peaks_per_image)
        return avg

    def peaks_per_spot(self):
        peaks_per_spot = []
        for i in self.images.keys():
            peaks_per_spot += self.images[i][1].peaks_per_spot(image=i)
        return peaks_per_spot

    def avg_peaks_per_spot(self):
        pps = self.peaks_per_spot()
        return sum(pps) / len(pps)

    def groups_per_image(self):
        groups_per_image = []
        for i in self.images.keys():
            groups_per_image.append(self.images[i][1].spot_count)
        return groups_per_image

    def avg_groups_per_image(self):
        groups = self.groups_per_image()
        return sum(groups) / len(groups)

    @classmethod
    def create_from_files(cls, inpath):
        images: Dict[Image, ImageDictValue] = {}
        if isinstance(inpath, str):
            inpath = glob.glob(inpath)
        for im in inpath:
            reader = Reader(im)
            image = Image(reader)
            imgstats = ImgStats(image)
            peakpositions = PeakPositions(image)
            images[image] = (peakpositions, imgstats)
        return cls(images)


class ImageDictStats:
    pass
