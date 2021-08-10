#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""
from peak_stats.reader.peaks import Image, Group, Peak
from scipy.spatial import ConvexHull
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Histograms global options
Figsize = (30, 20)
LabelFontsize = 35
TitleFontsize = 40
TickFontsize = 30
Facecolor = "g"
Alpha = 0.75
Dpi = 300

# 3D plots global options
Color = "red"


class ImgStats:

    def __init__(self, image: Image):
        self.group_count = image.group_count()
        self.peak_count = image.peak_count()
        self.avg_group_sigma_xyz = self.add_groups_sigma(image)
        self.avg_peaks = self.avg_peaks_per_group(image)
        self.photons = self.photons_in_image(image=image)
        self.photons_per_group = self.average_photons_per_group(image)

    @staticmethod
    def add_groups_sigma(image: Image):
        """Calculate sigma for x y and z"""
        avg_groups_sigma_z = []
        avg_groups_sigma_x = []
        avg_groups_sigma_y = []
        for group in image.groups:
            new_group_stat = GroupStats(group)
            avg_groups_sigma_x.append(new_group_stat.group_avg_sigma_x_pos_full)
            avg_groups_sigma_y.append(new_group_stat.group_avg_sigma_y_pos_full)
            avg_groups_sigma_z.append(new_group_stat.group_avg_sigma_z)
        return [avg_groups_sigma_x, avg_groups_sigma_y, avg_groups_sigma_z]

    @staticmethod
    def add_peak_sigma(image: Image):
        """Return peaks sigma in format[[simgax_peka1,...],[simgay_peak1...],[sigmaz_peak1]]"""
        sigma_x = []
        sigma_y = []
        sigma_z = []
        for group in image.groups:
            group_stats = GroupStats(group=group)
            sigma_x += group_stats.list_peak_sigma_x(group=group)
            sigma_y += group_stats.list_peak_sigma_y(group=group)
            sigma_z += group_stats.list_peak_sigma_z(group=group)
        return [sigma_x, sigma_y, sigma_z]

    @staticmethod
    def peaks_per_group(image: Image):
        peaks_in_spots = []
        for spot in image.groups:
            peaks_in_spots.append(len(spot))
        return peaks_in_spots

    def avg_peaks_per_group(self, image: Image):
        peaks = self.peaks_per_group(image)
        avg = sum(peaks) / len(peaks)
        return avg

    def average_sigma_x(self):
        average_sigma_x = sum(self.avg_group_sigma_xyz[0]) / len(self.avg_group_sigma_xyz[0])
        return average_sigma_x

    def average_sigma_y(self):
        average_sigma_y = sum(self.avg_group_sigma_xyz[1]) / len(self.avg_group_sigma_xyz[1])
        return average_sigma_y

    def average_sigma_z(self):
        average_sigma_z = sum(self.avg_group_sigma_xyz[2]) / len(self.avg_group_sigma_xyz[2])
        return average_sigma_z

    def average_sigma(self):
        average = (sum(self.avg_group_sigma_xyz[2]) / len(self.avg_group_sigma_xyz[2]) + sum(
            self.avg_group_sigma_xyz[1]) / len(self.avg_group_sigma_xyz[1]) + sum(self.avg_group_sigma_xyz[0]) / len(
            self.avg_group_sigma_xyz[0])) / 3
        return average

    def plot_sigma_x_hist(self, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        plt.hist(self.avg_group_sigma_xyz[0], 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Groups', fontsize=LabelFontsize)
        plt.xlabel('Sigma X', fontsize=LabelFontsize)
        plt.title('Sigma X', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if save:
            outfile = outdir + "_sigma_X"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_sigma_y_hist(self, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        plt.hist(self.avg_group_sigma_xyz[1], 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Groups', fontsize=LabelFontsize)
        plt.xlabel('Sigma Y', fontsize=LabelFontsize)
        plt.title('Sigma Y', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if save:
            outfile = outdir + "_sigma_Y"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_sigma_z_hist(self, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        plt.hist(self.avg_group_sigma_xyz[2], 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Groups', fontsize=LabelFontsize)
        plt.xlabel('Sigma Z', fontsize=LabelFontsize)
        plt.title('Sigma Z', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if save:
            outfile = outdir + "_sigma_Z"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_peak_sigma_x(self, x_sigma, save=False, outdir=None, draw_sigma=15):
        plt.figure(figsize=Figsize)
        plt.hist(x_sigma, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Peaks', fontsize=LabelFontsize)
        plt.xlabel('Sigma X', fontsize=LabelFontsize)
        plt.title('Sigma X', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if draw_sigma:
            plt.axvline(x=draw_sigma, color="red")
        if save:
            outfile = outdir + "_sigma_X_peaks_g"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_peak_sigma_y_(self, y_sigma, save=False, outdir=None, draw_sigma=15):
        plt.figure(figsize=Figsize)
        plt.hist(y_sigma, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Peaks', fontsize=LabelFontsize)
        plt.xlabel('Sigma Y', fontsize=LabelFontsize)
        plt.title('Sigma Y', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if draw_sigma:
            plt.axvline(x=draw_sigma, color="red")
        if save:
            outfile = outdir + "_sigma_Y_peaks_g"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_peak_sigma_z(self, sigma_z, save=False, outdir=None, draw_sigma=15):
        plt.figure(figsize=Figsize)
        plt.hist(sigma_z, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Peaks', fontsize=LabelFontsize)
        plt.xlabel('Sigma Z', fontsize=LabelFontsize)
        plt.title('Sigma Z', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if draw_sigma:
            plt.axvline(x=draw_sigma, color="red")
        if save:
            outfile = outdir + "_sigma_Z_peaks_g"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_average_peak_sigma(self, sigma_z, sigma_y, sigma_x, save=False, outdir=None, draw_sigma=15):
        avg = [0] * len(sigma_x)
        for i in range(len(sigma_x)):
            avg[i] += sigma_x[i]
            avg[i] += sigma_y[i]
            avg[i] += sigma_z[i]
            avg[i] /= 3
        plt.figure(figsize=Figsize)
        plt.hist(avg, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Peaks', fontsize=70)
        plt.xlabel('Average sigma', fontsize=70)
        plt.title('Average Peak Sigma', fontsize=100)
        plt.tick_params(axis='both', which='major', labelsize=50)
        plt.grid(True)
        if draw_sigma:
            plt.axvline(x=draw_sigma, color="red")
        if save:
            outfile = outdir + "_average_peaks_g_v2"
            plt.savefig(outfile, dpi=100, format="png")
        else:
            plt.show()
        plt.close()

    def plot_average_sigma(self, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        data = np.array(self.avg_group_sigma_xyz)
        avg = np.average(data, axis=0)
        plt.hist(avg, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Groups', fontsize=LabelFontsize)
        plt.xlabel('Average Sigma', fontsize=LabelFontsize)
        plt.title('Groups Average Sigma', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if save:
            outfile = outdir + "_sigma_avg"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_peaks_per_group_histogram(self, image: Image, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        data = self.peaks_per_group(image=image)
        plt.hist(data, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Groups', fontsize=LabelFontsize)
        plt.xlabel('Number of Peaks', fontsize=LabelFontsize)
        plt.title('Peaks per Group histogram', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.minorticks_on()
        plt.grid(True)
        if save:
            outfile = outdir + "_peaks_per_group"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()

    def plot_photons_per_group(self, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        plt.hist(self.photons_per_group, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Groups', fontsize=LabelFontsize)
        plt.xlabel('Average Number of Photons', fontsize=LabelFontsize)
        plt.title('Average number of photons per peak in group', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if save:
            outfile = outdir + "_photons_per_group"
            plt.savefig(outfile, dpi=Dpi, formay="png")
        else:
            plt.show()
        plt.close()

    def plot_photons_per_peaks(self, save=False, outdir=None):
        plt.figure(figsize=Figsize)
        plt.hist(self.photons, 50, facecolor=Facecolor, alpha=Alpha)
        plt.ylabel('Number of Peaks', fontsize=LabelFontsize)
        plt.xlabel('Number of Photons', fontsize=LabelFontsize)
        plt.title('Number of photons per peak', fontsize=TitleFontsize)
        plt.tick_params(axis='both', which='major', labelsize=TickFontsize)
        plt.grid(True)
        if save:
            outfile = outdir + "_photons_per_image"
            plt.savefig(outfile, dpi=Dpi, format="png")
        else:
            plt.show()
        plt.close()
        pass

    def save_statistics(self, output, image: Image, sigma_threshold):
        with open(output, 'w') as out:
            peaks = PeakPositions(image=image, sigma_threshold=sigma_threshold)
            groups = GroupPeakStats(image=image)
            peaks.plot_convex_hull(show=False)
            out.write("Group Count:  {}".format(self.group_count) + "\n")
            out.write("Peak Count:   {}".format(self.peak_count) + "\n")
            out.write("Average peaks per group   {}".format(self.peak_count / self.group_count) + "\n")
            out.write("Average sigma X:  {} nm".format(self.average_sigma_x()) + "\n")
            out.write("Average sigma Y:  {} nm".format(self.average_sigma_y()) + "\n")
            out.write("Average sigma Z:  {} nm".format(self.average_sigma_z()) + "\n")
            out.write("Convex hull volume:   {}".format(peaks.hull_volume) + "\n")
            out.write("Convex hull area:    {}".format(peaks.hull_area) + "\n")
            out.write("Average number of photons per peak:  {}".format(sum(self.photons) / len(self.photons)) + "\n")
            out.write("Average number of photons per group peak: {}".format(groups.average_photons()) + "\n")
            out.write("Average sigma X per group peaks: {} nm".format(groups.averege_sigma_x()) + "\n")
            out.write("Average sigma Y per group peaks: {} nm".format(groups.averege_sigma_y()) + "\n")
            out.write("Average sigma Z per group peaks: {} nm".format(groups.averege_sigma_z()))

    def print_statistics(self, image: Image, sigma_threshold=None):
        peaks = PeakPositions(image=image, sigma_threshold=sigma_threshold)
        groups = GroupPeakStats(image=image)
        peaks.plot_convex_hull(show=False)
        print("Group Count:  {}".format(self.group_count))
        print("Peak Count:   {}".format(self.peak_count))
        print("Average peaks per group   {}".format(self.peak_count / self.group_count))
        print("Average sigma X:  {} nm".format(self.average_sigma_x()))
        print("Average sigma Y:  {} nm".format(self.average_sigma_y()))
        print("Average sigma Z:  {} nm".format(self.average_sigma_z()))
        print("Convex hull volume:   {}".format(peaks.hull_volume))
        print("Convex hull area:    {}".format(peaks.hull_area))
        print("Average number of photons per peak:  {}".format(sum(self.photons) / len(self.photons)))
        print("Average number of photons per group peak: {}".format(groups.average_photons()))
        print("Average sigma X per group peaks: {} nm".format(groups.averege_sigma_x()))
        print("Average sigma Y per group peaks: {} nm".format(groups.averege_sigma_y()))
        print("Average sigma Z per group peaks: {} nm".format(groups.averege_sigma_z()))

    @staticmethod
    def photons_in_image(image: Image):
        """Return a list of number of phontons in each peak in the image"""
        photons_per_peak = []
        for spot in image.groups:
            photons_per_peak += GroupStats.photons(group=spot)
        return photons_per_peak

    @staticmethod
    def average_photons_per_group(image: Image):
        """Return a list of average number of photons per spot"""
        photons_per_spot = []
        for spot in image.groups:
            avg_photons = sum(GroupStats.photons(group=spot)) / len(spot)
            photons_per_spot.append(avg_photons)
        return photons_per_spot


class GroupStats:

    def __init__(self, group: Group):
        self.peak_count = len(group.peaks)
        self.group_avg_sigma_x_pos_full = self.calculate_avg_sigma_x(group)
        self.group_avg_sigma_y_pos_full = self.calculate_avg_sigma_y(group)
        self.group_avg_sigma_z = self.calculate_avg_sigma_z(group)

    def list_peak_sigma_x(self, group: Group, pixel_size=133):
        sigma_x = []
        for i in group.peaks:
            sigma_x.append(float(i.data["Sigma X Pos Full"]) * pixel_size)
        return sigma_x

    def list_peak_sigma_y(self, group: Group, pixel_size=133):
        sigma_y = []
        for i in group.peaks:
            sigma_y.append(float(i.data["Sigma Y Pos Full"]) * pixel_size)
        return sigma_y

    def list_peak_sigma_z(self, group: Group):
        sigma_z = []
        for i in group.peaks:
            sigma_z.append(float(i.data["Sigma Z"]))
        return sigma_z

    @staticmethod
    def calculate_avg_sigma_x(group: Group, pixel_size=133):
        sum_sigma = 0
        for i in group.peaks:
            sum_sigma += (float(i.data["Sigma X Pos Full"]) * pixel_size)
        avg_sigma = sum_sigma / len(group.peaks)
        return avg_sigma

    @staticmethod
    def calculate_avg_sigma_y(group: Group, pixel_size=133):
        sum_sigma = 0
        for i in group.peaks:
            sum_sigma += (float(i.data["Sigma Y Pos Full"]) * pixel_size)
        avg_sigma = sum_sigma / len(group.peaks)
        return avg_sigma

    @staticmethod
    def calculate_avg_sigma_z(group: Group):
        sum_sigma = 0
        for i in group.peaks:
            sum_sigma += float(i.data["Sigma Z"])
        avg_sigma = sum_sigma / len(group.peaks)
        return avg_sigma

    @staticmethod
    def photons(group: Group):
        """"return list of number of photons per peak in given spot"""
        photons = []
        for peak in group.peaks:
            photon = peak.data["6 N Photons"]
            photons.append(photon)
        return photons


class PeakPositions:

    def __init__(self, image: Image, sigma_threshold=1000, minimize=False):
        self.peaks_positions = self.image_peak_positions(image, sigma_threshold)
        self.parse_peak_positions_to_nm()
        if minimize:
            self.minimize_xy()
        self.hull_volume = None
        self.hull_area = None

    @staticmethod
    def single_peak_position(peak: Peak, sigma_threshold, pixel_size=133):
        peak_x = peak.data["X Position"]
        peak_y = peak.data["Y Position"]
        peak_z = peak.data["Unwrapped Z"]
        peak_position = [peak_x, peak_y, peak_z]
        if sigma_threshold:
            if peak.data["Sigma X Pos Full"] * pixel_size < sigma_threshold and peak.data[
                "Sigma Y Pos Full"] * pixel_size < sigma_threshold and peak.data["Sigma Z"] < sigma_threshold:
                return peak_position
            else:
                return None
        return peak_position

    def group_peak_positions(self, group: Group, sigma_threshold):
        spot_peaks_positions = []
        for peak in group.peaks:
            peak_positions = self.single_peak_position(peak, sigma_threshold)
            if peak_positions:
                spot_peaks_positions.append(peak_positions)
        return spot_peaks_positions

    def image_peak_positions(self, image: Image, sigma_threshold):
        image_peaks_positions = []
        for spot in image.groups:
            spot_peaks_positions = self.group_peak_positions(spot, sigma_threshold)
            image_peaks_positions += spot_peaks_positions
        return np.array(image_peaks_positions)

    def parse_peak_positions_to_nm(self, pixel_size=133):
        """X/Y positions in ASCII file are in pixels and Z position is in nm. """
        self.peaks_positions[:, 0] *= pixel_size
        self.peaks_positions[:, 1] *= pixel_size

    def minimize_xy(self):
        minimal = np.amin(self.peaks_positions, axis=0)
        self.peaks_positions[:, 0] = self.peaks_positions[:, 0] - minimal[0]
        self.peaks_positions[:, 1] = self.peaks_positions[:, 1] - minimal[1]
        return minimal

    def plot_peak_positions(self, title=None, outpath=None):
        """Draw 3D plot of peak positions."""
        fig = plt.figure(figsize=(10, 7))
        ax = Axes3D(fig)
        ax.scatter(self.peaks_positions[:, 0], self.peaks_positions[:, 1], self.peaks_positions[:, 2], s=20,
                   color=Color)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel("X [nm]", fontsize=12, labelpad=5)
        ax.set_ylabel("Y [nm]", fontsize=12, labelpad=5)
        ax.set_zlabel("Z [nm]", fontsize=12)
        if title:
            plt.title(title, fontsize=15)
        if outpath:
            plt.savefig(outpath, dpi=300, format="png")
        else:
            plt.show()

    def plot_convex_hull(self, title=None, outpath=None, show=True):
        """Plot convex hull on peak positions """
        hull = ConvexHull(self.peaks_positions)
        fig = plt.figure(figsize=(10, 7))
        ax = Axes3D(fig)
        ax.scatter(self.peaks_positions[:, 0], self.peaks_positions[:, 1], self.peaks_positions[:, 2], s=20,
                   color=Color)
        ax.tick_params(axis='both', which='major', labelsize=12)
        for s in hull.simplices:
            s = np.append(s, s[0])
            ax.plot(self.peaks_positions[s, 0], self.peaks_positions[s, 1], self.peaks_positions[s, 2],
                    color="grey", linewidth=1)
        if title:
            plt.title(title, fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel("X [nm]", fontsize=12, labelpad=5)
        ax.set_ylabel("Y [nm]", fontsize=12, labelpad=5)
        ax.set_zlabel("Z [nm]", fontsize=12)
        if outpath:
            plt.savefig(outpath, dpi=300)
        elif outpath is None and show is True:
            plt.show()
        plt.close()
        print(hull.volume)
        print(hull.area)
        self.hull_volume = hull.volume
        self.hull_area = hull.area


class GroupPeakStats:

    def __init__(self, image: Image):
        self.group_peaks = self.add_group_peaks(image=image)
        self.group_number = len(self.group_peaks)
        self.photons_per_group = self.photons_per_group_peak()
        self.positions = self.group_peaks_positions()
        self.minimize_group_peaks()

    @staticmethod
    def add_group_peaks(image: Image):
        group_peaks = []
        for spot in image.groups:
            group_peaks.append(spot.group_peak)
        return group_peaks

    # -----------------------------PHOTONS--------------------------------------
    def photons_per_group_peak(self):
        """Return a list of number of photons for each group peak."""
        photons = []
        for group_peak in self.group_peaks:
            photons.append(group_peak.photon_number)
        return photons

    def average_photons(self):
        """Return average number of photons per group"""
        return sum(self.photons_per_group) / len(self.photons_per_group)

    # ------------------------------XYZ POSITIONS---------------------------------
    def group_peaks_positions(self):
        """Return a list of group peaks xyz positions."""
        peak_positions = []
        for group_peak in self.group_peaks:
            position = [group_peak.x_position, group_peak.y_position, group_peak.z_position]
            peak_positions.append(position)
        return np.array(peak_positions)

    # -------------------------------SIGMA---------------------------------------

    def groups_sigma_x(self):
        """Return a list with x sigmas for all group peaks"""
        sigmas = []
        for group_peak in self.group_peaks:
            sigmas.append(group_peak.group_sigma_x)
        return sigmas

    def groups_sigma_y(self):
        """Return a list with y sigmas for all group peaks"""
        sigmas = []
        for group_peak in self.group_peaks:
            sigmas.append(group_peak.group_sigma_y)
        return sigmas

    def groups_sigma_z(self):
        """Return a list with z sigmas for all group peaks"""
        sigmas = []
        for group_peak in self.group_peaks:
            sigmas.append(group_peak.group_sigma_z)
        return sigmas

    def averege_sigma_x(self):
        """Return average group X sigma for the image"""
        sigmas = self.groups_sigma_x()
        return sum(sigmas) / len(sigmas)

    def averege_sigma_y(self):
        """Return average group Y sigma for the image"""
        sigmas = self.groups_sigma_y()
        return sum(sigmas) / len(sigmas)

    def averege_sigma_z(self):
        """Return average group Z sigma for the image"""
        sigmas = self.groups_sigma_z()
        return sum(sigmas) / len(sigmas)

    def minimize_group_peaks(self):
        """Move the structure to 0, 0 """
        minimal = np.amin(self.positions, axis=0)
        self.positions[:, 0] = self.positions[:, 0] - minimal[0]
        self.positions[:, 1] = self.positions[:, 1] - minimal[1]
        return minimal

    # -----------------PLOTTING----------------------------------

    def plot_3d_group_peaks(self, title="Group Peaks", outpath=None):
        """Draw 3D plot of group peaks in XYZ space."""
        fig = plt.figure(figsize=(10, 7))
        ax = Axes3D(fig)
        ax.scatter(self.positions[:, 0], self.positions[:, 1], self.positions[:, 2], s=20,
                   color=Color)
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_xlabel("X [nm]", fontsize=20, labelpad=10)
        ax.set_ylabel("Y [nm]", fontsize=20, labelpad=10)
        ax.set_zlabel("Z [nm]", fontsize=20, labelpad=10)
        if title:
            plt.title(title, fontsize=40)
        if outpath:
            plt.savefig(outpath, dpi=100, format="png")
        else:
            plt.show()

    def plot_group_peaks_convex_hull(self, title="Group Peaks Convex Hull", outpath=None, show=True):
        """Draw 3D plot of convex hull od group peaks and return the volume"""
        hull = ConvexHull(self.positions)
        fig = plt.figure(figsize=(10, 7))
        ax = Axes3D(fig)
        ax.scatter(self.positions[:, 0], self.positions[:, 1], self.positions[:, 2], s=20,
                   color=Color)
        ax.tick_params(axis='both', which='major', labelsize=12)
        for s in hull.simplices:
            s = np.append(s, s[0])
            ax.plot(self.positions[s, 0], self.positions[s, 1], self.positions[s, 2],
                    color="grey", linewidth=1)
        if title:
            plt.title(title, fontsize=40)
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_xlabel("X [nm]", fontsize=20, labelpad=10)
        ax.set_ylabel("Y [nm]", fontsize=20, labelpad=10)
        ax.set_zlabel("Z [nm]", fontsize=20, labelpad=10)
        if outpath:
            plt.savefig(outpath, dpi=100)
        elif outpath is None and show is True:
            plt.show()
        plt.close()
        print(hull.volume)
        print(hull.area)
