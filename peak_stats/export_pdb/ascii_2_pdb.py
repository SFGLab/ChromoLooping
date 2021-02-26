#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
from peak_stats.reader.peaks import Image
from peak_stats.statistics.single_img_stats import PeakPositions
import numpy as np


class Points:

    def __init__(self, image: Image, sigma_threshold):
        self.image = image
        self.points = self.all_peaks(sigma_threshold=sigma_threshold)
        self.minimal = self.minimize_xy()

    @staticmethod
    def save_pdb(filename, points):
        """Save pdb file with points"""
        atoms = ''
        n = len(points)
        for i in range(n):
            x = points[i][0]
            y = points[i][1]
            try:
                z = points[i][2]
            except IndexError:
                z = 0.0
            atoms += (
                '{0:6}{1:>5}  {2:3}{3:}{4:3} {5:}{6:>4}{7:}   {8:8.2f}{9:8.2f}{10:8.3f}{11:6.2f}{12:6.2f}{13:>12}\n'.format(
                    "ATOM", i + 1, 'B', ' ', 'BEA', 'A', i + 1, ' ', x, y, z, 0, 0, 'C'))
        pdb_file_content = atoms
        open(filename, 'w').write(pdb_file_content)

    def find_centroid_in_group(self, peaks_in_spots):
        """Given a list of lists with spot positions: [[peak1, peak2...],...[peak1, peak2,...]]
         return a list of central peaks for all spots"""
        spot_centroids = []
        for spot in peaks_in_spots:
            cetr = spot.mean(axis=0)
            spot_centroids.append(cetr)
        return spot_centroids

    def parse_image(self, sigma_threshold):
        """Find centroids for all groups in image."""
        groups_peaks = self.group_peaks(sigma_threshold=sigma_threshold)
        centroids = self.find_centroid_in_group(peaks_in_spots=groups_peaks)
        return centroids

    def all_peaks(self, sigma_threshold):
        """Create a list of all peaks."""
        peaksp = PeakPositions(self.image, sigma_threshold=sigma_threshold)
        return peaksp.peaks_positions

    def group_peaks(self, sigma_threshold, pixel_size=133):
        """Apply pixel size. Default X,Y values from PeakSelector are in voxels and for Z are in nm. Also move the image
         to the beginning of the coordinate system """
        peaks_per_spots = []
        peksp = PeakPositions(self.image)
        for spot in self.image.groups:
            positions = np.array(peksp.group_peak_positions(group=spot, sigma_threshold=sigma_threshold))
            if len(positions) > 0:
                positions[:, 0] *= pixel_size
                positions[:, 1] *= pixel_size
                positions[:, 0] = positions[:, 0] - self.minimal[0]
                positions[:, 1] = positions[:, 1] - self.minimal[1]
                peaks_per_spots.append(positions)
        return peaks_per_spots

    def minimize_xy(self):
        """Move the image to x=0 and y=0."""
        minimal = np.amin(self.points, axis=0)
        self.points[:, 0] = self.points[:, 0] - minimal[0]
        self.points[:, 1] = self.points[:, 1] - minimal[1]
        return minimal
