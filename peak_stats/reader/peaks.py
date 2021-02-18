#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
from peak_stats.reader.read_peaks import Reader
"""
In ASCII files from PeakSelector Z position is saved in nanometers, but X and Y positions are saved in pixels and have 
to be multiplied by pixel size.
"""

PIXEL_SIZE = 133

# todo - add some description to functions

class Peak:

    def __init__(self, data_row, header):
        self.data = {}
        for i in range(len(header)):
            self.data[header[i]] = float(data_row[i])


class Spot:

    def __len__(self):
        return len(self.peaks)

    def __init__(self, group_index, group_peak):
        self.group_index = group_index
        self.peaks = []
        self.group_peak = group_peak

    def add_peak(self, peak: Peak):
        self.peaks.append(peak)


class Image:

    def __init__(self, reader: Reader):
        self.fields = reader.head
        self.spots = []
        self.file_path = reader.path
        counter = 0
        for row in reader.lines:
            counter += 1
            self.add_peak(row, reader.head)

    def __hash__(self):
        return hash(self.file_path)

    def __eq__(self, other):
        return isinstance(other, Image) and self.file_path == other.file_path

    def spot_count(self):
        return len(self.spots)

    def peak_count(self):
        count = 0
        for spot in self.spots:
            count += len(spot)
        return count

    def add_spot(self, spot):
        self.spots.append(spot)

    def attributes(self):
        return self.fields

    def add_peak(self, data_row, head):
        # create peak
        new_peak = Peak(data_row=data_row, header=head)
        for i in self.spots:
            # check if exists spot for this peak
            if i.group_index == new_peak.data['18 Grouped Index']:
                i.add_peak(peak=new_peak)
                return
                # if not create a new spot and add to image
        group_peak = GroupPeak(peak=new_peak)
        new_spot = Spot(group_index=new_peak.data['18 Grouped Index'], group_peak=group_peak)
        new_spot.add_peak(peak=new_peak)
        self.add_spot(spot=new_spot)
        return


class GroupPeak:
    def __init__(self, peak: Peak):
        self.x_position = peak.data["Group X Position"] * PIXEL_SIZE
        self.y_position = peak.data["Group Y Position"] * PIXEL_SIZE
        self.z_position = peak.data["Group Z Position"]
        self.group_sigma_x = peak.data["Group Sigma X Pos"] * PIXEL_SIZE
        self.group_sigma_y = peak.data["Group Sigma Y Pos"] * PIXEL_SIZE
        self.group_sigma_z = peak.data["Group Sigma Z"]
        self.photon_number = peak.data["Group N Photons"]
