# from __future__ import print_function, absolute_import, division

import logging
import multiprocessing
import traceback
import lcbuilder.eleanor
import sys

sys.modules['eleanor'] = sys.modules['lcbuilder.eleanor']
import eleanor
from lcbuilder.eleanor.targetdata import TargetData
from lcbuilder.eleanor_manager import EleanorManager
import warnings
from itertools import chain
import PIL
import batman
import scipy
import foldedleastsquares
import lcbuilder
import lightkurve
from lightkurve.periodogram import BoxLeastSquaresPeriodogram
from lightkurve import MPLSTYLE
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import wotan
import yaml
from astropy.timeseries.periodograms import BoxLeastSquares
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from lcbuilder.constants import CUTOUT_SIZE, LIGHTKURVE_CACHE_DIR, ELEANOR_CACHE_DIR
from lcbuilder.helper import LcbuilderHelper
from lcbuilder.lcbuilder_class import LcBuilder
from lcbuilder.objectinfo.MissionObjectInfo import MissionObjectInfo
from lcbuilder.objectinfo.preparer.MissionLightcurveBuilder import MissionLightcurveBuilder
from lcbuilder.photometry.aperture_extractor import ApertureExtractor
from lcbuilder.star.EpicStarCatalog import EpicStarCatalog
from lcbuilder.star.HabitabilityCalculator import HabitabilityCalculator
from lcbuilder.star.KicStarCatalog import KicStarCatalog
from lcbuilder.star.TicStarCatalog import TicStarCatalog
from lightkurve import TessLightCurve, TessTargetPixelFile, KeplerTargetPixelFile
from matplotlib.colorbar import Colorbar
from matplotlib import patches
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table
from astropy.io import ascii
import astropy.visualization as stretching
from scipy import stats, ndimage

from watson import constants
import watson.tpfplotterSub.tpfplotter as tpfplotter
import pandas as pd
import os
from math import ceil, floor
from copy import deepcopy

from watson.neighbours import CreateStarCsvInput, create_star_csv, NeighbourInput, get_neighbour_lc
from watson.report import Report


class Watson:
    """
    Provides transiting candidate vetting information like centroids and spaceship motion, momentum dumps, neighbours
    curves inspection and more to give a deeper insight on the quality of the candidate signal.
    """
    def __init__(self, object_dir):
        self.object_dir = os.getcwd() if object_dir is None else object_dir
        self.data_dir = self.object_dir + "/"

    def vetting(self, id, period, t0, duration, depth, ffi, sectors, rp_rstar=None, a_rstar=None, cpus=None,
                cadence=None, lc_file=None, lc_data_file=None, tpfs_dir=None, apertures_file=None,
                create_fov_plots=False, cadence_fov=None, ra=None, dec=None, transits_list=None,
                v=None, j=None, h=None, k=None, clean=True, transits_mask=None):
        """
        Launches the whole vetting procedure that ends up with a validation report
        :param id: the target star id
        :param period: the period of the candidate in days
        :param t0: the epoch in days
        :param duration: the duration of the transit of the candidate in minutes
        :param depth: the depth of the transit of the candidate in ppts
        :param ffi: flag to specify whether the data to be used is FFI or short cadence curve
        :param sectors: sectors/quarters/campaigns to be used
        :param rp_rstar: Rp / Rstar
        :param a_rstar: Semi-major axis / Rstar
        :param cpus: number of cpus to be used
        :param cadence: the cadence to be used to download data, in seconds
        :param lc_file: the file containing the curve
        :param lc_data_file: the file containing the raw curve and the motion, centroids and quality flags
        :param tpfs_dir: the directory containing the tpf files
        :param apertures_file: the file containing the map of sectors->apertures
        :param create_fov_plots: whether to generate Field Of View plots.
        :param cadence_fov: the cadence to use to download fov_plots
        :param ra: the RA to use to download fov_plots
        :param dec: the DEC to use to download fov_plots
        :param transits_list: a list of dictionaries with shape: {'t0': value, 'depth': value, 'depth_err': value}
        :param v: star V magnitude
        :param j: star J magnitude
        :param h: star H magnitude
        :param k: star K magnitude
        :param clean: whether to clean all the pngs created for the final pdfs
        :param transits_mask: array with shape [{P:period, T0:t0, D:d}, ...] to use for transits masking before vetting
        """
        logging.info("------------------")
        logging.info("Candidate info")
        logging.info("------------------")
        logging.info("Id: %.s", id)
        logging.info("Period (d): %.2f", period)
        logging.info("Epoch (d): %.2f", t0)
        logging.info("Duration (min): %.2f", duration)
        logging.info("Depth (ppt): %.2f", depth)
        logging.info("Rp_Rstar: %.4f", rp_rstar)
        logging.info("a_Rstar: %.2f", a_rstar)
        logging.info("FFI: %s", ffi)
        logging.info("Sectors: %s", sectors)
        lc_builder = LcBuilder()
        if transits_mask is None:
            transits_mask = []
        if rp_rstar is None:
            rp_rstar = np.sqrt(depth / 1000)
        lc_build = None
        if lc_file is None or lc_data_file is None:
            lc_build = lc_builder.build(MissionObjectInfo(sectors, id, cadence=cadence,
                                                          initial_transit_mask=transits_mask), self.data_dir)
            lc_build.lc_data.to_csv(self.data_dir + "/lc_data.csv")
            lc_file = self.data_dir + "/lc.csv"
            lc_data_file = self.data_dir + "/lc_data.csv"
        if a_rstar is None and lc_build is None:
            raise ValueError("You need to define a_rstar if you are providing the lc_file and lc_data_file")
        if a_rstar is None:
            a_rstar = HabitabilityCalculator().calculate_semi_major_axis(period, lc_build.star_info.mass)
        if tpfs_dir is None:
            tpfs_dir = self.data_dir + "/tpfs/"
        if apertures_file is None:
            apertures_file = self.data_dir + "/apertures.yaml"
        index = 0
        vetting_dir = self.data_dir + "/vetting_" + str(index)
        while os.path.exists(vetting_dir) or os.path.isdir(vetting_dir):
            vetting_dir = self.data_dir + "/vetting_" + str(index)
            index = index + 1
        os.mkdir(vetting_dir)
        self.data_dir = vetting_dir
        try:
            transits_list_t0s, summary_list_t0s_indexes = self.__process(id, period, t0, duration, depth, rp_rstar, a_rstar,
                                                                 cpus, lc_file, lc_data_file, tpfs_dir,
                                                                 apertures_file, create_fov_plots, cadence_fov, ra,
                                                                 dec, transits_list, transits_mask)
            self.report(id, ra, dec, t0, period, duration, depth, transits_list_t0s, summary_list_t0s_indexes,
                        v, j, h, k, os.path.exists(tpfs_dir))
            if clean:
                for filename in os.listdir(self.data_dir):
                    if not filename.endswith(".pdf"):
                        os.remove(self.data_dir + "/" + filename)

        except Exception as e:
            traceback.print_exc()

    def report(self, id, ra, dec, t0, period, duration, depth, transits_list, summary_list_t0s_indexes, v, j, h, k,
               with_tpfs=True):
        file_name = "transits_validation_report.pdf"
        logging.info("Creating complete report")
        report = Report(self.data_dir, file_name, id, ra, dec, t0, period, duration, depth, transits_list, None,
                        v, j, h, k, with_tpfs)
        report.create_report()
        file_name = "transits_validation_report_summary.pdf"
        logging.info("Creating summary report")
        report = Report(self.data_dir, file_name, id, ra, dec, t0, period, duration, depth, transits_list,
                        summary_list_t0s_indexes, v, j, h, k, with_tpfs)
        report.create_report()


    def vetting_with_data(self, candidate_df, star, transits_df, cpus, create_fov_plots=False, cadence_fov=None,
                          transits_mask=None):
        """
        Same than vetting but receiving a candidate dataframe and a star dataframe with one row each.
        :param candidate_df: the candidate dataframe containing id, period, t0, transits and sectors data.
        :param star: the star dataframe with the star info.
        :param transits_df: a dataframe containing the transits information with columns 't0', 'depth' and 'depth_err'
        :param cpus: the number of cpus to be used.
        :param create_fov_plots: whether to generate Field Of View plots.
        :param cadence_fov: the cadence to use to download fov_plots
        :param transits_mask: array with shape [{P:period, T0:t0, D:d}, ...] to use for transits masking before vetting
        """
        if transits_mask is None:
            transits_mask = []
        df = candidate_df.iloc[0]
        # TODO get the transit time list
        id = df['id']
        period = df['period']
        t0 = df['t0']
        rp_rstar = df['rp_rs']
        a_rstar = df['a'] / star["R_star"] * constants.AU_TO_RSUN
        duration = df['duration']
        depth = df['depth']
        ffi = df['ffi']
        run = int(df['number'])
        curve = int(df['curve'])
        sectors = df['sectors']
        if isinstance(sectors, (int, np.integer)):
            sectors = [sectors]
        elif isinstance(sectors, (str)):
            sectors = sectors.split(',')
        lc_file = "/lc_" + str(curve) + ".csv"
        lc_file = self.object_dir + lc_file
        lc_data_file = self.object_dir + "/lc_data.csv"
        tpfs_dir = self.object_dir + "/tpfs"
        apertures_file = self.object_dir + "/apertures.yaml"
        try:
            self.vetting(id, period, t0, duration, depth, ffi, sectors, rp_rstar=rp_rstar, a_rstar=a_rstar, cpus=cpus,
                         lc_file=lc_file, lc_data_file=lc_data_file, tpfs_dir=tpfs_dir, apertures_file=apertures_file,
                         create_fov_plots=create_fov_plots, cadence_fov=cadence_fov, ra=star["ra"],
                         dec=star["dec"], transits_list=None if transits_df is None else transits_df.to_dict("list"),
                         transits_mask=transits_mask)
        except Exception as e:
            traceback.print_exc()

    def __process(self, id, period, t0, duration, depth, rp_rstar, a_rstar, cpus, lc_file, lc_data_file, tpfs_dir,
                  apertures_file, create_fov_plots=False, cadence_fov=None, ra_fov=None, dec_fov=None,
                  transits_list=None, transits_mask=None):
        """
        Performs the analysis to generate PNGs and Transits Validation Report.
        :param id: the target star id
        :param period: the period of the candidate in days
        :param t0: the epoch in days
        :param duration: the duration of the transit of the candidate in minutes
        :param depth: the depth of the transit of the candidate in ppts
        :param sectors: sectors/quarters/campaigns to be used
        :param rp_rstar: Rp / Rstar
        :param a_rstar: Semi-major axis / Rstar
        :param cpus: number of cpus to be used
        :param lc_file: the file containing the curve
        :param lc_data_file: the file containing the raw curve and the motion, centroids and quality flags
        :param tpfs_dir: the directory containing the tpf files
        :param apertures_file: the file containing the apertures
        :param create_fov_plots: whether to create FOV plots
        :param cadence_fov: the cadence to use to download fov_plots
        :param ra_fov: the RA to use to download fov_plots
        :param dec_fov: the DEC to use to download fov_plots
        :param transits_list: a list of dictionaries with shape: {'t0': value, 'depth': value, 'depth_err': value}
        :param transits_mask: array with shape [{P:period, T0:t0, D:d}, ...] to use for transits masking before vetting
        """
        logging.info("Running Transit Plots")
        lc, lc_data, lc_data_norm, tpfs = Watson.initialize_lc_and_tpfs(id, lc_file, lc_data_file, tpfs_dir,
                                                          transits_mask=transits_mask)
        apertures = None
        sectors = None
        if os.path.exists(apertures_file):
            apertures = yaml.load(open(apertures_file), yaml.SafeLoader)
            apertures = apertures["sectors"]
            sectors = [sector for sector in apertures.keys()]
            mission, mission_prefix, mission_int_id = LcBuilder().parse_object_info(id)
            if create_fov_plots:
                if cadence_fov is None:
                    cadence_fov = LcbuilderHelper.compute_cadence(lc.time.value)
                Watson.vetting_field_of_view(self.data_dir, mission, mission_int_id, cadence_fov, ra_fov, dec_fov,
                                             list(apertures.keys()), "tpf", apertures, cpus)
        summary_t0s_indexes = None
        if transits_list is not None:
            transits_list_not_nan_indexes = \
                Watson.plot_transits_statistics(self.data_dir, id, t0, period, transits_list)
            transit_t0s_list = np.array(transits_list["t0"])[transits_list_not_nan_indexes]
            transit_depths = np.array(transits_list["depth"])[transits_list_not_nan_indexes]
            summary_t0s_indexes = np.argwhere((transit_depths == np.max(transit_depths)) |
                                              (transit_depths == np.min(transit_depths))).flatten()
            if len(transit_depths) > 2:
                closest_depths_to_mean = np.abs(transit_depths - depth)
                summary_t0s_indexes = np.append(summary_t0s_indexes, np.argmin(closest_depths_to_mean))
        else:
            last_time = lc.time.value[len(lc.time.value) - 1]
            first_time = lc.time.value[0]
            num_of_transits_back = int(floor(((t0 - first_time) / period)))
            transits_lists_back = t0 - period * np.arange(num_of_transits_back, 0, -1) if num_of_transits_back > 0 else np.array([])
            num_of_transits = int(ceil(((last_time - t0) / period)))
            transit_lists = t0 + period * np.arange(0, num_of_transits)
            transit_lists = np.append(transits_lists_back, transit_lists)
            time_as_array = lc.time.value
            plot_range = duration / 3600 * 2
            transits_in_data = [time_as_array[(transit > time_as_array - plot_range) & (transit < time_as_array + plot_range)] for
                                transit in transit_lists]
            transit_t0s_list = transit_lists[[len(transits_in_data_set) > 0 for transits_in_data_set in transits_in_data]]
        mission, mission_prefix, target_id = MissionLightcurveBuilder().parse_object_id(id)
        if (ra_fov is not None and dec_fov is not None):
            Watson.plot_folded_tpfs(self.data_dir, mission_prefix, mission, target_id, ra_fov, dec_fov, lc, lc_data,
                                    tpfs, lc_file, lc_data_file, tpfs_dir, sectors, period, t0, duration, depth / 1000,
                                    rp_rstar, a_rstar, transits_mask, transit_t0s_list, cpus)
        self.plot_folded_curve(self.data_dir, id, lc, period, t0, duration, depth / 1000, rp_rstar, a_rstar)
        Watson.plot_all_folded_cadences(self.data_dir, mission_prefix, mission, target_id, lc, sectors, period, t0,
                                        duration, depth / 1000, rp_rstar, a_rstar, cpus)
        #self.plot_nb_stars(self.data_dir, mission, id, lc, period, t0, duration, depth / 1000, cpus)
        plot_transits_inputs = []
        for index, transit_times in enumerate(transit_t0s_list):
            plot_transits_inputs.append(SingleTransitProcessInput(self.data_dir, str(id), index, lc_file, lc_data_file,
                                                                  tpfs_dir, apertures, transit_times, depth / 1000,
                                                                  duration, period, rp_rstar, a_rstar, transits_mask))
        with multiprocessing.Pool(processes=cpus) as pool:
            pool.map(Watson.plot_single_transit, plot_transits_inputs)
        return transit_t0s_list, summary_t0s_indexes

    @staticmethod
    def initialize_lc_and_tpfs(id, lc_file, lc_data_file, tpfs_dir, transits_mask=None):
        if transits_mask is None:
            transits_mask = []
        mission, mission_prefix, mission_int_id = LcBuilder().parse_object_info(id)
        lc = pd.read_csv(lc_file, header=0)
        lc_data = None
        if os.path.exists(lc_data_file):
            lc_data = pd.read_csv(lc_data_file, header=0)
            lc_data_norm = Watson.normalize_lc_data(lc_data)
        time, flux, flux_err = lc["#time"].values, lc["flux"].values, lc["flux_err"].values
        for transit_mask in transits_mask:
            logging.info('* Transit mask with P=%.2f d, T0=%.2f d, Dur=%.2f min *', transit_mask["P"],
                         transit_mask["T0"], transit_mask["D"])
            mask = foldedleastsquares.transit_mask(time, transit_mask["P"], transit_mask["D"] / 60 / 24,
                                                   transit_mask["T0"])
            time = time[~mask]
            flux = flux[~mask]
            flux_err = flux_err[~mask]
        lc = TessLightCurve(time=time, flux=flux, flux_err=flux_err, quality=np.zeros(len(time)))
        lc.extra_columns = []
        tpfs = []
        if os.path.exists(tpfs_dir):
            for tpf_file in sorted(os.listdir(tpfs_dir)):
                tpf = TessTargetPixelFile(tpfs_dir + "/" + tpf_file) if mission == lcbuilder.constants.MISSION_TESS else \
                    KeplerTargetPixelFile(tpfs_dir + "/" + tpf_file)
                for transit_mask in transits_mask:
                    mask = foldedleastsquares.transit_mask(tpf.time.value, transit_mask["P"], transit_mask["D"] / 60 / 24,
                                                           transit_mask["T0"])
                    tpf = tpf[~mask]
                tpfs.append(tpf)
        return lc, lc_data, lc_data_norm, tpfs

    @staticmethod
    def normalize_lc_data(lc_data):
        logging.info("Normalizing lc_data")
        lc_data_copy = lc_data.copy()
        time = lc_data_copy["time"].to_numpy()
        dif = time[1:] - time[:-1]
        jumps = np.where(np.abs(dif) > 0.2)[0]
        jumps = np.append(jumps, len(lc_data_copy))
        previous_jump_index = 0
        for jumpIndex in jumps:
            token = lc_data_copy["centroids_x"][previous_jump_index:jumpIndex]
            lc_data_copy["centroids_x"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            token = lc_data_copy["centroids_y"][previous_jump_index:jumpIndex]
            lc_data_copy["centroids_y"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            token = lc_data_copy["motion_x"][previous_jump_index:jumpIndex]
            lc_data_copy["motion_x"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            token = lc_data_copy["motion_y"][previous_jump_index:jumpIndex]
            lc_data_copy["motion_y"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            previous_jump_index = jumpIndex
        return lc_data_copy

    @staticmethod
    def plot_transits_statistics(data_dir, id, epoch, period, transits_list):
        fig, axs = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
        fig.suptitle(str(id) + ' Transits depth analysis T0=' + str(round(epoch, 2)) + ' P=' + str(round(period, 2)) + 'd', size=18)
        transits_list_not_nan_t0s_indexes = np.argwhere(~np.isnan(transits_list["t0"])).flatten()
        transits_list_not_nan_depths_indexes = np.argwhere(~np.isnan(transits_list["depth"])).flatten()
        transits_list_not_nan_indexes = np.intersect1d(transits_list_not_nan_t0s_indexes, transits_list_not_nan_depths_indexes)
        transits_list_t0s = np.array(transits_list["t0"])[transits_list_not_nan_indexes]
        transits_list_depths = np.array(transits_list["depth"])[transits_list_not_nan_indexes]
        transits_list_depths_err = np.array(transits_list["depth_err"])[transits_list_not_nan_indexes]
        even_transits_indexes = np.argwhere((np.abs((transits_list_t0s - epoch) % (2 * period)) < 0.1) |
                                            (np.abs((transits_list_t0s - epoch) % (2 * period)) > (
                                                        2 * period) - 0.1)).flatten()
        odd_transits_indexes = np.argwhere((np.abs((transits_list_t0s - epoch) % (2 * period)) > period - 0.05) &
                                           (np.abs((transits_list_t0s - epoch) % (2 * period)) < period + 0.05)).flatten()
        axs.axhline(y=np.mean(transits_list_depths), color='purple', alpha=0.3,
                    ls='-', lw=2, label='Depth Mean')
        axs.axhline(y=np.percentile(transits_list_depths, 84), color='purple', alpha=0.3,
                    ls='--', lw=2, label='Depth 1-sigma confidence')
        axs.axhline(y=np.percentile(transits_list_depths, 18), color='purple', alpha=0.3,
                    ls='--', lw=2)
        axs.axhline(y=np.mean(transits_list_depths[even_transits_indexes]), color='blue', alpha=0.3,
                    ls='-', lw=2, label='Depth Mean Even')
        axs.axhline(y=np.mean(transits_list_depths[odd_transits_indexes]), color='red', alpha=0.3,
                    ls='-', lw=2, label='Depth Mean Odd')
        axs.errorbar(x=even_transits_indexes, y=transits_list_depths[even_transits_indexes],
                     yerr=transits_list_depths_err[even_transits_indexes],
                     fmt="o", color="blue", ecolor="cyan", label="Even transits")
        axs.errorbar(x=odd_transits_indexes, y=transits_list_depths[odd_transits_indexes],
                     yerr=transits_list_depths_err[odd_transits_indexes],
                     fmt="o", color="red", ecolor="darkorange", label="Odd transits")
        axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axs.set_xlabel("Transit number")
        axs.set_ylabel("Depth (ppt)")
        plt.savefig(data_dir + "/transit_depths.png")
        plt.clf()
        plt.close()
        return transits_list_not_nan_indexes

    @staticmethod
    def plot_single_transit(single_transit_process_input):
        """
        Plots the single transit info: single transit focused curve, drift and background plots, small vs used aperture
        photometry and tpf flux values around the transit.
        @param single_transit_process_input: wrapper class to provide pickable inputs for multiprocessing
        """
        lc, lc_data, lc_data_norm, tpfs = Watson.initialize_lc_and_tpfs(single_transit_process_input.id,
                                                          single_transit_process_input.lc_file,
                                                          single_transit_process_input.lc_data_file,
                                                          single_transit_process_input.tpfs_dir,
                                                          transits_mask=single_transit_process_input.transits_mask)
        transit_time = single_transit_process_input.transit_times
        duration = single_transit_process_input.duration / 60 / 24
        fig = plt.figure(figsize=(24, 6), constrained_layout=True)
        fig.suptitle('Vetting of ' + str(single_transit_process_input.id) + ' single transit no. ' +
                     str(single_transit_process_input.index) +
                     ' at T0=' + str(round(transit_time, 2)) + 'd', size=26)
        gs = gridspec.GridSpec(2, 3, hspace=0.4, wspace=0.1)  # 2 rows, 3 columns
        ax1 = fig.add_subplot(gs[0, 0])  # First row, first column
        ax2 = fig.add_subplot(gs[0, 1])  # First row, second column
        ax3 = fig.add_subplot(gs[0, 2])  # First row, third column
        ax4 = fig.add_subplot(gs[1, 0])  # First row, third column
        ax5 = fig.add_subplot(gs[1, 1])  # First row, third column
        ax6 = fig.add_subplot(gs[1, 2])  # First row, third column
        axs = [ax1, ax2, ax3, ax4, ax5, ax6]
        sort_args = np.argsort(lc.time.value)
        time = lc.time.value[sort_args]
        flux = lc.flux.value[sort_args]
        flux_err = lc.flux_err.value[sort_args]
        plot_range = duration * 2
        zoom_mask = np.where((time > transit_time - plot_range) & (time < transit_time + plot_range))
        plot_time = time[zoom_mask]
        plot_flux = flux[zoom_mask]
        zoom_lc_data = None
        if lc_data_norm is not None:
            zoom_lc_data = lc_data_norm[(lc_data_norm["time"] > transit_time - plot_range) &
                                        (lc_data_norm["time"] < transit_time + plot_range)]
        aperture_mask = None
        chosen_aperture_lc = None
        smaller_aperture_lc = None
        eroded_aperture_mask = None
        tpf_short_framed = None
        if len(plot_time) > 0:
            for tpf in tpfs:
                tpf_zoom_mask = tpf[(tpf.time.value > transit_time - plot_range) & (tpf.time.value < transit_time +
                                                                                    plot_range)]
                if len(tpf_zoom_mask) > 0:
                    if hasattr(tpf, 'sector') and tpf.sector is not None:
                        sector = tpf.sector
                    elif hasattr(tpf, 'campaign') and tpf.campaign:
                        sector = tpf.campaign
                    else:
                        sector = tpf.quarter
                    tpf_short_framed = tpf[(tpf.time.value > transit_time - plot_range) &
                                           (tpf.time.value < transit_time + plot_range)]
                    if len(tpf_short_framed) == 0:
                        break
                    if not isinstance(single_transit_process_input.apertures, (dict)) and \
                            np.isnan(single_transit_process_input.apertures):
                        chosen_aperture_lc = lc
                        smaller_aperture_lc = lc
                        aperture_mask = [[]]
                    else:
                        aperture_mask = ApertureExtractor.from_pixels_to_boolean_mask(
                            single_transit_process_input.apertures[sector], tpf.column, tpf.row, tpf.shape[2], tpf.shape[1])
                        eroded_aperture_mask = ndimage.binary_erosion(aperture_mask)
                        chosen_aperture_lc = tpf.to_lightcurve(aperture_mask=aperture_mask)
                    if True in eroded_aperture_mask:
                        smaller_aperture_lc = tpf.to_lightcurve(aperture_mask=eroded_aperture_mask)
                    break
            single_transit_file = single_transit_process_input.data_dir + "/single_transit_" + \
                                  str(single_transit_process_input.index) + "_T0_" + str(transit_time) + ".png"
            tpf_single_transit_file = single_transit_process_input.data_dir + "/tpf_single_transit_" + \
                                      str(single_transit_process_input.index) + "_T0_" + str(transit_time) + ".png"
            t1 = transit_time - duration / 2
            t4 = transit_time + duration / 2
            model_time, model_flux = Watson.get_transit_model(duration, transit_time,
                                                              (transit_time - plot_range,
                                                               transit_time + plot_range),
                                                              single_transit_process_input.depth,
                                                              single_transit_process_input.period,
                                                              single_transit_process_input.rp_rstar,
                                                              single_transit_process_input.a_rstar)
            momentum_dumps_lc_data = None
            if zoom_lc_data is not None and not zoom_lc_data["quality"].isnull().all():
                momentum_dumps_lc_data = zoom_lc_data[np.bitwise_and(zoom_lc_data["quality"].to_numpy(),
                                                                     constants.MOMENTUM_DUMP_QUALITY_FLAG) >= 1]

            axs[0].plot(model_time, model_flux, color="red")
            axs[0].scatter(plot_time, plot_flux, color="darkorange", label="Photometry used aperture")
            axs[0].set_xlim([transit_time - plot_range, transit_time + plot_range])
            axs[0].set_title("Single Transit")
            if momentum_dumps_lc_data is not None and len(momentum_dumps_lc_data) > 0:
                first_momentum_dump = True
                for index, momentum_dump_row in momentum_dumps_lc_data.iterrows():
                    if first_momentum_dump:
                        axs[0].axvline(x=momentum_dump_row["time"], color='purple', alpha=0.3,
                                       ls='--', lw=2, label='Momentum dump')
                        first_momentum_dump = False
                    else:
                        axs[0].axvline(x=momentum_dump_row["time"], color='purple', alpha=0.3,
                                       ls='--', lw=2)
            axs[0].legend(loc='upper left')
            axs[0].set_xlim([transit_time - plot_range, transit_time + plot_range])
            axs[0].set_xlabel("Time (d)")
            axs[0].set_ylabel("Flux norm.")
            if zoom_lc_data is not None:
                axs[1].scatter(zoom_lc_data["time"], zoom_lc_data["motion_x"], color="red",
                               label="X-axis motion")
                axs[1].scatter(zoom_lc_data["time"], zoom_lc_data["centroids_x"],
                               color="black", label="X-axis centroids")
                axs[1].axvline(x=transit_time - duration / 2, color='r', label='T1')
                axs[1].axvline(x=transit_time + duration / 2, color='r', label='T4')
                axs[1].legend(loc='upper left')
            axs[1].set_xlim([transit_time - plot_range, transit_time + plot_range])
            axs[1].set_title("X-axis drift")
            axs[1].set_xlabel("Time (d)")
            axs[1].set_ylabel("Normalized Y-axis data")
            if zoom_lc_data is not None:
                axs[2].scatter(zoom_lc_data["time"], zoom_lc_data["motion_y"], color="red",
                               label="Y-axis motion")
                axs[2].scatter(zoom_lc_data["time"], zoom_lc_data["centroids_y"],
                               color="black", label="Y-axis centroids")
                axs[2].axvline(x=t1, color='r', label='T1')
                axs[2].axvline(x=t4, color='r', label='T4')
                axs[2].legend(loc='upper left')
            axs[2].set_xlim([transit_time - plot_range, transit_time + plot_range])
            axs[2].set_title("Y-axis drift")
            axs[2].set_xlabel("Time (d)")
            axs[2].set_ylabel("Normalized Y-axis data")
            if smaller_aperture_lc is not None:
                axs[3].plot(model_time, model_flux, color="red")
                axs[3].set_xlim([transit_time - plot_range, transit_time + plot_range])
                axs[3].set_title("SAP comparison")
                axs[3].set_xlim([transit_time - plot_range, transit_time + plot_range])
                axs[3].set_xlabel("Time (d)")
                axs[3].set_ylabel("Flux norm.")
                chosen_aperture_lc.flux = wotan.flatten(chosen_aperture_lc.time.value,
                                                        chosen_aperture_lc.flux.value, window_length=0.75,
                                                        return_trend=False, method="biweight", break_tolerance=0.5)
                smaller_aperture_lc.flux = wotan.flatten(smaller_aperture_lc.time.value,
                                                         smaller_aperture_lc.flux.value, window_length=0.75,
                                                         return_trend=False, method="biweight", break_tolerance=0.5)
                chosen_aperture_lc = chosen_aperture_lc[
                    (chosen_aperture_lc.time.value - plot_range < transit_time) &
                    (chosen_aperture_lc.time.value + plot_range > transit_time)]
                smaller_aperture_lc = smaller_aperture_lc[
                    (smaller_aperture_lc.time.value - plot_range < transit_time) &
                    (smaller_aperture_lc.time.value + plot_range > transit_time)]
                axs[3].scatter(chosen_aperture_lc.time.value, chosen_aperture_lc.flux.value, color="darkorange",
                               label="Photometry used aperture")
                axs[3].scatter(smaller_aperture_lc.time.value, smaller_aperture_lc.flux.value,
                               color="c", label="Photometry smaller aperture")
            axs[3].legend(loc='upper left')
            if tpf_short_framed is not None:
                axs[4] = tpf_short_framed.plot(axs[4], aperture_mask=aperture_mask)
                axs[4].set_title("TPF apertures comparison")
                if smaller_aperture_lc is not None:
                    parsed_aperture = tpf_short_framed._parse_aperture_mask(eroded_aperture_mask)
                    for i in range(tpf_short_framed.shape[1]):
                        for j in range(tpf_short_framed.shape[2]):
                            if parsed_aperture[i, j]:
                                rect = patches.Rectangle(
                                    xy=(j + tpf_short_framed.column - 0.5, i + tpf_short_framed.row - 0.5),
                                    width=1,
                                    height=1,
                                    color='black',
                                    fill=False,
                                    hatch="\\\\",
                                )
                                axs[4].add_patch(rect)
            if zoom_lc_data is not None:
                axs[5].scatter(zoom_lc_data["time"], zoom_lc_data["background_flux"],
                               color="blue", label="Background Flux")
                axs[5].axvline(x=t1, color='r', label='T1')
                axs[5].axvline(x=t4, color='r', label='T4')
                axs[5].legend(loc='upper left')
            axs[5].set_xlim([transit_time - plot_range, transit_time + plot_range])
            axs[5].set_title("Background flux")
            axs[5].set_xlabel("Time (d)")
            axs[5].set_ylabel("Background flux (e/s)")
            plt.savefig(single_transit_file, dpi=100, bbox_inches='tight')
            plt.clf()
            plt.close()
            if tpf_short_framed is not None:
                tpf_short_framed.plot_pixels(aperture_mask=aperture_mask, markersize=1)
                plt.savefig(tpf_single_transit_file, dpi=100)
                plt.clf()
                plt.close()
                images_list = [single_transit_file, tpf_single_transit_file]
                imgs = [PIL.Image.open(i) for i in images_list]
                imgs[0] = imgs[0].resize((imgs[1].size[0],
                                          int(imgs[1].size[0] / imgs[0].size[0] * imgs[0].size[1])),
                                          PIL.Image.ANTIALIAS)
                img_merge = np.vstack((np.asarray(i) for i in imgs))
                img_merge = PIL.Image.fromarray(img_merge)
                img_merge.save(single_transit_file, quality=95, optimize=True)
                os.remove(tpf_single_transit_file)
            logging.info("Processed single transit plot for T0=%.2f", transit_time)
        else:
            logging.info("Not plotting single transit for T0=%.2f as the data is empty", transit_time)

    @staticmethod
    def plot_folded_curve(file_dir, id, lc, period, epoch, duration, depth, rp_rstar, a_rstar):
        """
        Plots the phase-folded curve of the candidate for period, 2 * period and period / 2.
        @param file_dir: the directory to store the plot
        @param id: the target id
        @param period: the transit period
        @param epoch: the transit epoch
        @param duration: the transit duration
        @param depth: the transit depth
        """
        duration = duration / 60 / 24
        figsize = (16, 16)  # x,y
        rows = 3
        cols = 2
        fig, axs = plt.subplots(rows, cols, figsize=figsize, constrained_layout=True)
        logging.info("Preparing folded light curves for target")
        #TODO bins = None for FFI
        bins = 100
        Watson.compute_phased_values_and_fill_plot(id, axs[0][0], lc, period, epoch, depth, duration, rp_rstar, a_rstar,
                                                   bins=bins)
        Watson.compute_phased_values_and_fill_plot(id, axs[0][1], lc, period, epoch + period / 2, depth, duration,
                                                   rp_rstar, a_rstar, bins=bins)
        period = 2 * period
        Watson.compute_phased_values_and_fill_plot(id, axs[1][0], lc, period, epoch, depth, duration, rp_rstar, a_rstar,
                                                   bins=bins)
        Watson.compute_phased_values_and_fill_plot(id, axs[1][1], lc, period, epoch + period / 2, depth, duration,
                                                   rp_rstar, a_rstar, bins=bins)
        period = period / 4
        Watson.compute_phased_values_and_fill_plot(id, axs[2][0], lc, period, epoch, depth, duration, rp_rstar, a_rstar,
                                                   bins=bins)
        Watson.compute_phased_values_and_fill_plot(id, axs[2][1], lc, period, epoch + period / 2, depth, duration,
                                                   rp_rstar, a_rstar, bins=bins)
        plt.savefig(file_dir + "/odd_even_folded_curves.png", dpi=200)
        fig.clf()
        plt.close(fig)

    @staticmethod
    def plot_all_folded_cadences(file_dir, mission_prefix, mission, id, lc, sectors, period, epoch, duration, depth, rp_rstar,
                                 a_rstar, cpus=os.cpu_count() - 1):
        updated_eleanor = False
        bins = 100
        fig, axs = plt.subplots(3, 1, figsize=(15, 10))
        duration = duration / 60 / 24
        duration_to_period = duration / period
        cadences = ['fast', 'short', 'long']
        for index, cadence in enumerate(cadences):
            lc = None
            found_sectors = []
            axs[index].set_title(mission_prefix + " " + str(id) + " " + str(found_sectors) + ": " + cadence)
            if mission == lcbuilder.constants.MISSION_TESS and cadence == 'long':
                author = "TESS-SPOC"
            elif mission == lcbuilder.constants.MISSION_TESS and cadence != 'long':
                author = "SPOC"
            elif mission == lcbuilder.constants.MISSION_KEPLER:
                author = "Kepler"
            elif mission == lcbuilder.constants.MISSION_K2:
                author = "K2"
            if mission == "TESS" and cadence == 'long':
                lcs = lightkurve.search_lightcurve(
                    mission_prefix + " " + str(id),
                    mission=mission,
                    sector=sectors,
                    campaign=sectors,
                    quarter=sectors,
                    author=author,
                    cadence=cadence
                ).download_all(download_dir=os.path.expanduser('~') + '/' + LIGHTKURVE_CACHE_DIR)
                if lcs is None:
                    if not updated_eleanor:
                        EleanorManager.update()
                        updated_eleanor = True
                    star = eleanor.multi_sectors(tic=round(id), sectors='all',
                                                 post_dir=os.path.expanduser('~') + '/' + ELEANOR_CACHE_DIR,
                                                 metadata_path=os.path.expanduser('~') + '/' + ELEANOR_CACHE_DIR)
                    data = []
                    for s in star:
                        datum = TargetData(s, height=CUTOUT_SIZE, width=CUTOUT_SIZE, do_pca=True)
                        data.append(datum)
                    lc = data[0].to_lightkurve(data[0].pca_flux).remove_nans().flatten()
                    if len(data) > 1:
                        for datum in data[1:]:
                            lc = lc.append(datum.to_lightkurve(datum.pca_flux).remove_nans().flatten())
                    lc = lc.remove_nans()
                else:
                    matching_objects = []
                    for i in range(0, len(lcs.data)):
                        if lcs.data[i].label == "TIC " + str(id):
                            if lc is None:
                                lc = lcs.data[i].normalize()
                            else:
                                lc = lc.append(lcs.data[i].normalize())
                        else:
                            matching_objects.append(lcs.data[i].label)
                    if lc is None:
                        continue
                    else:
                        if mission == lcbuilder.constants.MISSION_TESS:
                            found_sectors = lcs.sector
                        elif mission == lcbuilder.constants.MISSION_KEPLER:
                            found_sectors = lcs.quarter
                        elif mission == lcbuilder.constants.MISSION_K2:
                            found_sectors = lcs.campaign
                    lc = lc.remove_nans()
                    if mission == lcbuilder.constants.MISSION_K2:
                        lc = lc.to_corrector("sff").correct(windows=20)
            else:
                lcs = lightkurve.search_lightcurve(
                    mission_prefix + " " + str(id),
                    mission=mission,
                    sector=sectors,
                    campaign=sectors,
                    quarter=sectors,
                    author=author,
                    cadence=cadence
                ).download_all(download_dir=os.path.expanduser('~') + '/' + LIGHTKURVE_CACHE_DIR)
                if lcs is None:
                    continue
                matching_objects = []
                for i in range(0, len(lcs.data)):
                    if lcs.data[i].label == mission_prefix + " " + str(id):
                        if lc is None:
                            lc = lcs.data[i].normalize()
                        else:
                            lc = lc.append(lcs.data[i].normalize())
                    else:
                        matching_objects.append(lcs.data[i].label)
                if lc is None:
                    continue
                else:
                    if mission == lcbuilder.constants.MISSION_TESS:
                        found_sectors = lcs.sector
                    elif mission == lcbuilder.constants.MISSION_KEPLER:
                        found_sectors = lcs.quarter
                    elif mission == lcbuilder.constants.MISSION_K2:
                        found_sectors = lcs.campaign
                lc = lc.remove_nans()
                lc = lc.remove_outliers(sigma_lower=float('inf'), sigma_upper=3)
                if mission == lcbuilder.constants.MISSION_K2:
                    lc = lc.to_corrector("sff").correct(windows=20)
            Watson.compute_phased_values_and_fill_plot(id, axs[index], lc, period, epoch, depth, duration,
                                                       rp_rstar, a_rstar, bins=bins)
            axs[index].set_title(mission_prefix + " " + str(id) + " " + str(found_sectors) + ": " + cadence)
        file = file_dir + '/folded_cadences.png'
        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
        plt.savefig(file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        plt.clf()

    @staticmethod
    def plot_folded_tpf(fold_tpf_input):
        lc_file = fold_tpf_input['lc_file']
        lc_data_file = fold_tpf_input['lc_data_file']
        tpfs_dir = fold_tpf_input['tpfs_dir']
        transits_mask = fold_tpf_input['transits_mask']
        period = fold_tpf_input['period']
        duration = fold_tpf_input['duration']
        t0s_list = fold_tpf_input['t0s_list']
        epoch = fold_tpf_input['epoch']
        mission = fold_tpf_input['mission']
        mission_prefix = fold_tpf_input['mission_prefix']
        id = fold_tpf_input['id']
        file_dir = fold_tpf_input['file_dir']
        lc, lc_data, lc_data_norm, tpfs = Watson.initialize_lc_and_tpfs(mission_prefix + ' ' + str(id), lc_file,
                                                                        lc_data_file, tpfs_dir,
                                                                        transits_mask=transits_mask)
        tpf = tpfs[fold_tpf_input['index']]
        pixel_values_i = np.array(range(tpf[0].shape[1]))
        pixel_values_j = np.array(range(tpf[0].shape[2]))
        tpf_lc_data = lc_data[(lc_data['time'] >= tpf.time.value[0]) & (lc_data['time'] <= tpf.time.value[-1])].dropna()
        sector_name, sector = LcbuilderHelper.mission_lightkurve_sector_extraction(mission, tpf)
        logging.info("Computing TPF centroids for %s %.0f", sector_name, sector)
        cadence_s = np.round(np.nanmedian(tpf.time.value[1:] - tpf.time.value[:-1]) * 3600 * 24)
        cadences_per_transit = LcbuilderHelper.estimate_transit_cadences(cadence_s, duration * 2)
        t0s_in_tpf_indexes = np.argwhere((t0s_list > tpf.time.value[0] - duration) &
                                         (t0s_list < tpf.time.value[-1] + duration)).flatten()
        if len(t0s_in_tpf_indexes) == 0:
            logging.warning("No transit was present in %s %.0f", sector_name, sector)
            return None, None, None
        tpf_t0s_list = t0s_list[t0s_in_tpf_indexes]
        good_quality_t0s = []
        for t0 in tpf_t0s_list:
            t0_in_tpf_indexes = \
                np.argwhere((tpf.time.value > t0 - duration) & (tpf.time.value < t0 + duration)).flatten()
            cadences_ratio = len(t0_in_tpf_indexes) / cadences_per_transit
            if cadences_ratio >= 0.75:
                good_quality_t0s.append(t0)
            else:
                tpf = tpf[(tpf.time.value < t0 - duration) | (tpf.time.value > t0 + duration)]
                tpf_lc_data = tpf_lc_data[(tpf_lc_data['time'] < t0 - duration) | (tpf_lc_data['time'] > t0 + duration)]
        if len(good_quality_t0s) == 0:
            logging.warning("There were transits T0s in %s %.0f but they had no good quality", sector_name, sector)
            return None, None, None
        snr_map, ax = Watson.plot_pixels(tpf, title=mission_prefix + ' ' + str(id) + ' ' + sector_name + ' ' +
                                                    str(sector) + ' TPF BLS Analysis',
                                         period=period, epoch=epoch, duration=duration, aperture_mask="pipeline")
        plt.savefig(file_dir + '/folded_tpf_' + str(sector) + '.png', dpi=200, bbox_inches='tight')
        plt.clf()
        hdu = tpf.hdu[2].header
        wcs = WCS(hdu)
        centroids_coords = np.array([[coord.ra.value, coord.dec.value] for coord in
                                     wcs.pixel_to_world(tpf_lc_data['centroids_x'], tpf_lc_data['centroids_y'])])
        centroids_offsets_ra = wotan.flatten(tpf_lc_data['time'].to_numpy(), centroids_coords[:, 0], duration * 4)
        # we sum 180 to the declination because flatten does not tolerate negative values
        centroids_offsets_dec = wotan.flatten(tpf_lc_data['time'].to_numpy(), centroids_coords[:, 1] + 180, duration * 4)
        source_offset = Watson.light_centroid(snr_map, pixel_values_i, pixel_values_j)
        source_offset = wcs.pixel_to_world(source_offset[1], source_offset[0])
        time = foldedleastsquares.core.fold(tpf.time.value, period, epoch + period / 2)
        lc_df = pd.DataFrame(columns=['time', 'flux', 'time_folded'])
        lc_df['time'] = tpf.time.value
        lc_df['time_folded'] = time
        lc_df_it = lc_df.loc[(lc_df['time_folded'] >= 0.5 - duration / 2) &
                             (lc_df['time_folded'] <= 0.5 + duration / 2)]
        lc_df_oot = lc_df.loc[
            ((lc_df['time_folded'] < 0.5 - duration / 2) & (lc_df['time_folded'] > 0.5 - 3 * duration / 2)) |
            ((lc_df['time_folded'] > 0.5 + duration / 2) & (lc_df['time_folded'] < 0.5 + 3 * duration / 2))]
        tpf_sub = np.zeros((tpf.shape[1], tpf.shape[2]))
        for i in np.arange(0, tpf.shape[1]):
            for j in np.arange(0, tpf.shape[2]):
                pixel_flux = wotan.flatten(tpf.time.value, tpf.flux[:, i, j], duration * 4, method='biweight')
                lc_df = pd.DataFrame(columns=['time', 'flux', 'time_folded'])
                lc_df['time'] = tpf.time.value
                lc_df['time_folded'] = time
                lc_df['flux'] = pixel_flux
                lc_df = lc_df.sort_values(by=['time_folded'], ascending=True)
                lc_df_it = lc_df.loc[(lc_df['time_folded'] >= 0.5 - duration / 2) &
                                     (lc_df['time_folded'] <= 0.5 + duration / 2)]
                lc_df_oot = lc_df.loc[
                    ((lc_df['time_folded'] < 0.5 - duration / 2) & (lc_df['time_folded'] > 0.5 - 3 * duration / 2)) |
                    ((lc_df['time_folded'] > 0.5 + duration / 2) & (lc_df['time_folded'] < 0.5 + 3 * duration / 2))]
                tpf_fluxes_oot = lc_df_oot['flux'].to_numpy()
                tpf_fluxes_it = lc_df_it['flux'].to_numpy()
                tpf_sub[i, j] = (np.nanmedian(tpf_fluxes_oot) - np.nanmedian(tpf_fluxes_it)) / \
                                np.sqrt((np.nanstd(tpf_fluxes_oot) ** 2) + (np.nanstd(tpf_fluxes_it) ** 2))
        hdu = tpf.hdu[2].header
        wcs = WCS(hdu)
        light_centroid_sub = Watson.light_centroid(tpf_sub, pixel_values_i, pixel_values_j)
        light_centroids_sub_coord = wcs.pixel_to_world(light_centroid_sub[1], light_centroid_sub[0])
        return (source_offset.ra.value, source_offset.dec.value), \
               (light_centroids_sub_coord.ra.value, light_centroids_sub_coord.dec.value), \
               (tpf_lc_data['time'].to_numpy(), centroids_offsets_ra, centroids_offsets_dec)

    @staticmethod
    def plot_folded_tpfs(file_dir, mission_prefix, mission, id, ra, dec, lc, lc_data, tpfs, lc_file, lc_data_file,
                         tpfs_dir, sectors, period, epoch, duration, depth, rp_rstar, a_rstar, transits_mask, t0s_list,
                         cpus=os.cpu_count() - 1):
        duration = duration / 60 / 24
        duration_to_period = duration / period
        i0 = tpfs[0].shape[1] // 2
        j0 = tpfs[0].shape[2] // 2
        yvec = np.array(range(150))
        target_coords = (ra, dec)
        logging.info("Computing TPF centroids")
        source_offsets = []
        tpf_fold_inputs = []
        for index, tpf in enumerate(tpfs):
            tpf_fold_inputs.append({'id': id, 'tpfs_dir': tpfs_dir, 'lc_file': lc_file, 'lc_data_file': lc_data_file,
                                    'index': index, 'duration': duration, 't0s_list': t0s_list, 'period': period,
                                    'epoch': epoch, 'transits_mask': transits_mask, 'mission': mission,
                                    'mission_prefix': mission_prefix, 'file_dir': file_dir})
        with multiprocessing.Pool(processes=cpus) as pool:
            results_fold = pool.map(Watson.plot_folded_tpf, tpf_fold_inputs)
        light_centroids_sub = []
        light_centroids_sub_coords = []
        centroids_offsets_ra_list = []
        centroids_offsets_dec_list = []
        centroids_offsets_time_list = []
        for index, tpf in enumerate(tpfs):
            if results_fold[index][0] is not None:
                source_offsets.append((results_fold[index][0][0], results_fold[index][0][1]))
            if results_fold[index][1] is not None:
                light_centroids_sub_coords.append((results_fold[index][1][0], results_fold[index][1][1]))
            if results_fold[index][2] is not None:
                centroids_offsets_time_list.append(results_fold[index][2][0])
                centroids_offsets_ra_list.append(results_fold[index][2][1])
                centroids_offsets_dec_list.append(results_fold[index][2][2])
        # TODO we don't manage to get a nice plot from this
        # centroid_coords_df = pd.DataFrame(columns=['time', 'time_folded', 'centroids_ra', 'centroids_dec'])
        # centroid_coords_df['time'] = list(chain.from_iterable(centroids_offsets_time_list))
        # centroid_coords_df['centroids_ra'] = list(chain.from_iterable(centroids_offsets_ra_list))
        # centroid_coords_df['centroids_dec'] = list(chain.from_iterable(centroids_offsets_dec_list))
        # centroid_coords_df['time_folded'] = foldedleastsquares.fold(centroid_coords_df['time'].to_numpy(), period, epoch + period / 2)
        # centroid_coords_df = centroid_coords_df[(centroid_coords_df['time_folded'] > 0.5 - duration_to_period * 2) &
        #                                         (centroid_coords_df['time_folded'] < 0.5 + duration_to_period * 2)]
        # centroid_coords_df = centroid_coords_df.sort_values(by=['time_folded'], ascending=True)
        # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
        # axs[0].scatter(centroid_coords_df['time_folded'], centroid_coords_df['centroids_ra'])
        # axs[1].scatter(centroid_coords_df['time_folded'], centroid_coords_df['centroids_dec'])
        # fig.suptitle(mission_prefix + ' ' + str(id) + ' - Transit centroid offsets')
        # plt.show()
        light_centroids_sub_ra = np.nanmedian(np.array(light_centroids_sub_coords)[:, 0])
        light_centroids_sub_dec = np.nanmedian(np.array(light_centroids_sub_coords)[:, 1])
        light_centroids_sub_ra_err = np.nanstd(np.array(light_centroids_sub_coords)[:, 0])
        light_centroids_sub_dec_err = np.nanstd(np.array(light_centroids_sub_coords)[:, 1])
        source_offset_ra = np.nanmedian(np.array(source_offsets)[:, 0])
        source_offset_dec = np.nanmedian(np.array(source_offsets)[:, 1])
        source_offset_ra_err = np.nanstd(np.array(source_offsets)[:, 0])
        source_offset_dec_err = np.nanstd(np.array(source_offsets)[:, 1])
        offset_ra = np.mean([source_offset_ra, light_centroids_sub_ra])
        offset_dec = np.mean([source_offset_dec, light_centroids_sub_dec])
        offset_ra_err = np.sqrt(source_offset_ra_err ** 2 + light_centroids_sub_ra_err ** 2)
        offset_dec_err = np.sqrt(source_offset_dec_err ** 2 + light_centroids_sub_dec_err ** 2)
        tpf = tpfs[0]
        hdu = tpf.hdu[2].header
        wcs = WCS(hdu)
        offset_px = wcs.all_world2pix(offset_ra, offset_dec, 0)
        light_centroids_sub_offset_px = wcs.all_world2pix(light_centroids_sub_ra, light_centroids_sub_dec, 0)
        source_offset_px = wcs.all_world2pix(source_offset_ra, source_offset_dec, 0)
        c1 = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
        c2 = SkyCoord(ra=offset_ra * u.degree, dec=offset_dec * u.degree, frame='icrs')
        distance_sub_arcs = c1.separation(c2).value * 60 * 60
        target_pixels = wcs.all_world2pix(ra, dec, 0)
        ax = tpf.plot(aperture_mask='pipeline')
        ax.plot([tpf.column + target_pixels[0]], [tpf.row + target_pixels[1]], marker="*", markersize=14,
                color="blue", label='target')
        offset_err = offset_ra_err if offset_ra_err > offset_dec_err else offset_dec_err
        offset_err = offset_err * 60 * 60
        offset_px_err = offset_err / LcbuilderHelper.mission_pixel_size(mission)
        circle1 = plt.Circle((tpf.column + offset_px[0], tpf.row + offset_px[1]),
                             offset_px_err, color='orange', fill=False)
        ax.add_patch(circle1)
        ax.plot([tpf.column + offset_px[0]], [tpf.row + offset_px[1]], marker="o",
                markersize=10, color="red", label='Diff image offset')
        ax.plot([tpf.column + source_offset_px[0]], [tpf.row + source_offset_px[1]], marker="*",
                markersize=4, color="green", label='Diff image offset')
        ax.plot([tpf.column + light_centroids_sub_offset_px[0]], [tpf.row + light_centroids_sub_offset_px[1]],
                marker="*", markersize=4, color="cyan", label='Diff image offset')
        ax.set_title(mission_prefix + ' ' + str(id) + " Source offsets - " +
                     str(np.round(distance_sub_arcs, 2)) + r'$\pm$' + str(np.round(offset_err, 2)) + "''")
        plt.savefig(file_dir + '/source_offsets.png', dpi=200, bbox_inches='tight')
        return source_offset_ra, source_offset_dec, distance_sub_arcs

    @staticmethod
    def light_centroid(snr_map, pixel_values_i, pixel_values_j):
        snr_i_0 = 0
        snr_j_0 = 0
        for i in pixel_values_i:
            for j in pixel_values_j:
                snr_i_0 = snr_i_0 + (snr_map[i, j] ** 2) * i
                snr_j_0 = snr_j_0 + (snr_map[i, j] ** 2) * j
        snr_div = 0
        for i in pixel_values_i:
            for j in pixel_values_j:
                snr_div = snr_div + (snr_map[i, j] ** 2)
        c_i = snr_i_0 / snr_div
        c_j = snr_j_0 / snr_div
        #mass_center = ndimage.measurements.center_of_mass(snr_map)
        return c_i, c_j

    @staticmethod
    def plot_pixels(
            tpf,
            ax=None,
            periodogram=False,
            aperture_mask=None,
            show_flux=False,
            corrector_func=None,
            style="lightkurve",
            title=None,
            markersize=0.5,
            period=None,
            epoch=None,
            duration=1,
            **kwargs,
    ):
        if style == "lightkurve" or style is None:
            style = MPLSTYLE
        if corrector_func is None:
            corrector_func = lambda x: x.remove_outliers()
        if show_flux:
            cmap = plt.get_cmap()
            norm = plt.Normalize(
                vmin=np.nanmin(tpf.flux[0].value), vmax=np.nanmax(tpf.flux[0].value)
            )
        mask = tpf._parse_aperture_mask(aperture_mask)

        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore", category=(RuntimeWarning)
            )

            # get an aperture mask for each pixel
            masks = np.zeros(
                (tpf.shape[1] * tpf.shape[2], tpf.shape[1], tpf.shape[2]),
                dtype="bool",
            )
            for i in range(tpf.shape[1] * tpf.shape[2]):
                masks[i][np.unravel_index(i, (tpf.shape[1], tpf.shape[2]))] = True

            pixel_list = []
            pixel_model_list = []
            lc = None
            bls_results = np.zeros((tpf.shape[1], tpf.shape[2])).tolist()
            for j in range(tpf.shape[1] * tpf.shape[2]):
                lc = tpf.to_lightcurve(aperture_mask=masks[j])
                lc = corrector_func(lc)
                lc.flux = wotan.flatten(lc.time.value, lc.flux.value, duration * 4, method='biweight')
                if period is not None:
                    duration_to_period = duration / period
                    lc_df = pd.DataFrame(columns=['time', 'time_folded', 'flux', 'flux_err'])
                    lc_df['time'] = lc.time.value
                    lc_df['time_folded'] = foldedleastsquares.core.fold(lc.time.value, period, epoch + period / 2)
                    lc_df['flux'] = lc.flux.value
                    lc_df['flux_err'] = lc.flux_err.value
                    lc_df = lc_df.sort_values(by=['time_folded'], ascending=True)
                    lc_df = lc_df[(lc_df['time_folded'] > 0.5 - duration_to_period * 3) & (lc_df['time_folded'] < 0.5 + duration_to_period * 3)]
                    lc = TessLightCurve(time=lc_df['time_folded'], flux=lc_df['flux'], flux_err=lc_df['flux_err'])
                    bls = BoxLeastSquares(lc_df['time_folded'].to_numpy(), lc_df['flux'].to_numpy(), lc_df['flux_err'].to_numpy())
                    result = bls.power([1],
                                       np.linspace(duration_to_period - duration_to_period / 2, duration_to_period * 3 / 2, 10))
                    x, y = np.unravel_index(j, (tpf.shape[1], tpf.shape[2]))
                    bls_results[x][y] = result
                if periodogram:
                    try:
                        pixel_list.append(lc.to_periodogram(**kwargs))
                    except IndexError:
                        pixel_list.append(None)
                else:
                    if len(lc.remove_nans().flux) == 0:
                        pixel_list.append(None)
                        pixel_model_list.append(None)
                    else:
                        pixel_list.append(lc)
                        if period is not None:
                            model = np.ones(len(lc))
                            it_mask = np.argwhere((lc.time.value > 0.5 - duration_to_period / 2) & (lc.time.value < 0.5 + duration_to_period / 2)).flatten()
                            model[it_mask] = 1 - result['depth'][0]
                            pixel_model_list.append(model)
        with plt.style.context(style):
            if ax is None:
                fig = plt.figure()
                ax = plt.gca()
                set_size = True
            else:
                fig = ax.get_figure()
                set_size = False

            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            if periodogram:
                ax.set(
                    title=title,
                    xlabel="Frequency / Column (pixel)",
                    ylabel="Power / Row (pixel)",
                )
            else:
                ax.set(
                    title=title,
                    xlabel="Time / Column (pixel)",
                    ylabel="Flux / Row (pixel)",
                )

            gs = gridspec.GridSpec(
                tpf.shape[1], tpf.shape[2], wspace=0.01, hspace=0.01
            )

            for k in range(tpf.shape[1] * tpf.shape[2]):
                if pixel_list[k]:
                    x, y = np.unravel_index(k, (tpf.shape[1], tpf.shape[2]))

                    # Highlight aperture mask in red
                    if aperture_mask is not None and mask[x, y]:
                        rc = {"axes.linewidth": 4, "axes.edgecolor": "purple"}
                    else:
                        rc = {"axes.linewidth": 1}
                    with plt.rc_context(rc=rc):
                        gax = fig.add_subplot(gs[tpf.shape[1] - x - 1, y])

                    # Determine background and foreground color
                    if show_flux:
                        gax.set_facecolor(cmap(norm(tpf.flux.value[0, x, y])))
                        markercolor = "white"
                    else:
                        markercolor = "black"

                    # Plot flux or periodogram
                    if periodogram:
                        gax.plot(
                            pixel_list[k].frequency.value,
                            pixel_list[k].power.value,
                            marker="None",
                            color=markercolor,
                            lw=markersize,
                        )
                    else:
                        gax.plot(
                            pixel_list[k].time.value,
                            pixel_list[k].flux.value,
                            marker=".",
                            color=markercolor,
                            ms=markersize,
                            lw=0,
                        )
                        if period is not None:
                            gax.plot(
                                pixel_list[k].time.value,
                                pixel_model_list[k],
                                marker=".",
                                color='red',
                                alpha=0.8,
                                ms=markersize,
                                lw=0,
                            )

                    gax.margins(y=0.1, x=0)
                    gax.set_xticklabels("")
                    gax.set_yticklabels("")
                    gax.set_xticks([])
                    gax.set_yticks([])

                    # add row/column numbers to start / end
                    if x == 0 and y == 0:
                        gax.set_xlabel(f"{tpf.column}")
                        gax.set_ylabel(f"{tpf.row}")
                    if x == 0 and y == tpf.shape[2] - 1:  # lower right
                        gax.set_xlabel(f"{tpf.column + tpf.shape[2] - 1}")
                    if x == tpf.shape[1] - 1 and y == 0:  # upper left
                        gax.set_ylabel(f"{tpf.row + tpf.shape[1] - 1}")

            if set_size:  # use default size when caller does not supply ax
                fig.set_size_inches((y * 1.5, x * 1.5))
        transit_times_score = np.zeros((tpf.shape[1], tpf.shape[2])).tolist()
        duration_score = np.zeros((tpf.shape[1], tpf.shape[2])).tolist()
        depth_score = np.zeros((tpf.shape[1], tpf.shape[2])).tolist()
        residuals = np.zeros((tpf.shape[1], tpf.shape[2])).tolist()
        for k in range(tpf.shape[1] * tpf.shape[2]):
            if pixel_list[k] is None:
                residuals[x][y] = np.inf
                continue
            x, y = np.unravel_index(k, (tpf.shape[1], tpf.shape[2]))
            max_power_index = np.argwhere(bls_results[x][y].power == np.nanmax(bls_results[x][y].power)).flatten()[0]
            best_epoch = bls_results[x][y].transit_time[max_power_index]
            best_duration = bls_results[x][y].duration[max_power_index]
            best_power = bls_results[x][y].power[max_power_index]
            best_depth = bls_results[x][y].depth[max_power_index]
            residuals[x][y] = np.sqrt(np.sum((pixel_list[k].flux.value - pixel_model_list[k]) ** 2))
            transit_times_score[x][y] = best_epoch
            duration_score[x][y] = best_duration
            depth_score[x][y] = best_depth
        transit_times_score = np.array(transit_times_score)
        duration_score = np.array(duration_score)
        depth_score = np.array(depth_score)
        transit_times_score = 1 / np.abs(transit_times_score - 0.5)
        duration_score = 1 / np.abs(duration_score - duration_to_period)
        total_score = np.sqrt(transit_times_score * duration_score * depth_score / residuals)
        total_score = np.nan_to_num(total_score, nan=np.nanmedian(total_score))
        snr_map = total_score / np.std(total_score)
        return snr_map, ax

    @staticmethod
    def plot_nb_stars(file_dir, mission, id, lc, period, epoch, duration, depth, cores=os.cpu_count()):
        if mission == lcbuilder.constants.MISSION_TESS:
            pixel_size = 20.25
            star_catalog = TicStarCatalog()
            author = "TESS-SPOC"
        elif mission == lcbuilder.constants.MISSION_KEPLER:
            star_catalog = KicStarCatalog()
            pixel_size = 4
            author = "Kepler"
        elif mission == lcbuilder.constants.MISSION_K2:
            star_catalog = EpicStarCatalog()
            pixel_size = 4
            author = "K2"
        search_radius = lcbuilder.constants.CUTOUT_SIZE / 2
        star_csv_file = \
            create_star_csv(CreateStarCsvInput(None, mission, id, pixel_size, search_radius, None, star_catalog,
                                               file_dir))
        plot_grid_size = 4
        fig, axs = plt.subplots(plot_grid_size, plot_grid_size, figsize=(16, 16))
        stars_df = pd.read_csv(star_csv_file)
        stars_df = stars_df.loc[~np.isnan(stars_df['id'])]
        stars_df = stars_df.sort_values(by=['dist_arcsec'], ascending=True)
        page = 0
        file = file_dir + '/star_nb_' + str(page) + '.png'
        duration = duration / 60 / 60
        if mission == "TESS":
            EleanorManager.update()
        neighbour_inputs = []
        for index, star_row in stars_df.iterrows():
            neighbour_inputs.append(
                NeighbourInput(index, mission, author, round(star_row['id']), period, epoch, duration))
        with multiprocessing.Pool(processes=1) as pool:
            lc_dfs = pool.map(get_neighbour_lc, neighbour_inputs)
        for index, neighbour_input in enumerate(neighbour_inputs):
            star_dist = stars_df.loc[neighbour_input.index, 'dist_arcsec']
            lc_df = lc_dfs[index]
            if len(lc_df) == 0:
                continue
            axs[index // 4 % 4][index % 4].set_title(str(id) + " - " + str(np.round(star_dist, 2)))
            axs[index // 4 % 4][index % 4].scatter(lc_df['folded_time'], lc_df['flux'], color='black', alpha=0.5)
            bins = 20
            if len(lc_df) > bins:
                bin_means, bin_edges, binnumber = stats.binned_statistic(lc_df['folded_time'], lc_df['flux'],
                                                                         statistic='mean', bins=bins)
                bin_width = (bin_edges[1] - bin_edges[0])
                bin_centers = bin_edges[1:] - bin_width / 2
                axs[index // 4 % 4][index % 4].scatter(bin_centers, bin_means, color='orange')
            #axs[index // 4][index % 4].axhline(1 - depth, color="red")
            if index % (plot_grid_size * plot_grid_size) == 0 and index > 0:
                plt.savefig(file, dpi=200)
                plt.close(fig)
                plt.clf()
                file = file_dir + '/star_nb_' + str(page) + '.png'
                if index + 1 < len(stars_df):
                    fig, axs = plt.subplots(plot_grid_size, plot_grid_size, figsize=(16, 16))
            page = index // (plot_grid_size * plot_grid_size)
        if index % (plot_grid_size * plot_grid_size) != 0:
            plt.savefig(file, dpi=200)
            plt.close(fig)
            plt.clf()
        return

    @staticmethod
    def compute_phased_values_and_fill_plot(id, axs, lc, period, epoch, depth, duration, rp_rstar, a_rstar, range=5,
                                            bins=None, bin_err_mode="flux_err"):
        """
        Phase-folds the input light curve and plots it centered in the given epoch
        @param id: the candidate name
        @param axs: the plot axis to be drawn
        @param lc: the lightkurve object containing the data
        @param period: the period for the phase-folding
        @param epoch: the epoch to center the fold
        @param depth: the transit depth
        @param duration: the transit duration
        @param range: the range to be used from the midtransit time in half-duration units.
        @param bins: the number of bins
        @params bin_err_mode: either 'bin' or 'flux_err' for flux_err std.
        @return: the drawn axis and the computed bins
        """
        time = foldedleastsquares.core.fold(lc.time.value, period, epoch + period / 2)
        axs.scatter(time, lc.flux.value, 2, color="blue", alpha=0.1)
        sort_args = np.argsort(time)
        time = time[sort_args]
        flux = lc.flux.value[sort_args]
        flux_err = lc.flux_err.value[sort_args]
        half_duration_phase = duration / 2 / period
        folded_plot_range = half_duration_phase * range
        folded_plot_range = folded_plot_range if folded_plot_range < 0.5 else 0.5
        folded_phase_zoom_mask = np.where((time > 0.5 - folded_plot_range) &
                                          (time < 0.5 + folded_plot_range))
        folded_phase = time[folded_phase_zoom_mask]
        folded_y = flux[folded_phase_zoom_mask]
        folded_y_err = flux_err[folded_phase_zoom_mask]
        axs.set_xlim([0.5 - folded_plot_range, 0.5 + folded_plot_range])
        # TODO if FFI no binning
        if bins is not None and len(folded_y) > bins:
            bin_means, bin_edges, binnumber = stats.binned_statistic(folded_phase, folded_y,
                                                                     statistic='mean', bins=bins)
            bin_width = (bin_edges[1] - bin_edges[0])
            bin_centers = bin_edges[1:] - bin_width / 2
            if bin_err_mode == 'flux_err':
                bin_stds, _, _ = stats.binned_statistic(folded_phase, folded_y_err, statistic='mean', bins=bins)
            else:
                bin_stds, _, _ = stats.binned_statistic(folded_phase, folded_y, statistic='std', bins=bins)
            bin_nan_args = np.isnan(bin_stds)
            axs.errorbar(bin_centers[~bin_nan_args], bin_means[~bin_nan_args],
                         yerr=bin_stds[~bin_nan_args] / 2, xerr=bin_width / 2, marker='o', markersize=2,
                         color='darkorange', alpha=1, linestyle='none')
        else:
            bin_centers = folded_phase
            bin_means = folded_y
            bin_stds = folded_y_err
        model_time, model_flux = Watson.get_transit_model(half_duration_phase * 2, 0.5,
                                                          (0.5 - half_duration_phase * range, 0.5 + half_duration_phase * range),
                                                          depth, period, rp_rstar, a_rstar, 2 * len(time))
        axs.plot(model_time, model_flux, color="red")
        axs.set_title(str(id) + " Folded Curve with P={:.2f}d and T0={:.2f}".format(period, epoch))
        axs.set_xlabel("Time (d)")
        axs.set_ylabel("Flux norm.")
        if len(folded_y) > 0:
            axs.set_ylim(np.min(folded_y), np.max(folded_y))
        #axs.set_ylim([1 - 3 * depth, 1 + 3 * depth])
        logging.info("Processed phase-folded plot for P=%.2f and T0=%.2f", period, epoch)
        return axs, bin_centers, bin_means, bin_stds

    @staticmethod
    #TODO build model from selected transit_template
    def get_transit_model(duration, t0, start_end, depth, period, rp_to_rstar, a_to_rstar, model_len=10000):
        t = np.linspace(-6, 6, model_len)
        ma = batman.TransitParams()
        ma.t0 = 0  # time of inferior conjunction
        ma.per = 365  # orbital period, use Earth as a reference
        ma.rp = rp_to_rstar  # planet radius (in units of stellar radii)
        ma.a = a_to_rstar  # semi-major axis (in units of stellar radii)
        ma.inc = 90  # orbital inclination (in degrees)
        ma.ecc = 0  # eccentricity
        ma.w = 0  # longitude of periastron (in degrees)
        ma.u = [0.4804, 0.1867]  # limb darkening coefficients
        ma.limb_dark = "quadratic"  # limb darkening model
        m = batman.TransitModel(ma, t)  # initializes model
        model = m.light_curve(ma)  # calculates light curve
        model_intransit = np.argwhere(model < 1)[:, 0]
        model_time = np.linspace(start_end[0], start_end[1], len(model))
        in_transit_indexes = np.where((model_time > t0 - duration / 2) & (model_time < t0 + duration / 2))[0]
        model_time_in_transit = model_time[in_transit_indexes]
        scaled_intransit = np.interp(
            np.linspace(model_time_in_transit[0], model_time_in_transit[-1], len(in_transit_indexes)),
            model_time[model_intransit], model[model_intransit])
        model = np.full((model_len), 1.0)
        model[in_transit_indexes] = scaled_intransit
        model[model < 1] = 1 - ((1 - model[model < 1]) * depth / (1 - np.min(model)))
        return model_time, model

    @staticmethod
    def plot_tpf(tpf, sector, aperture, dir):
        logging.info("Plotting FOV curves for sector %.0f", sector)
        if not os.path.exists(dir):
            os.mkdir(dir)
        tpf.plot_pixels(aperture_mask=aperture)
        plt.savefig(dir + "/fov_Flux_pixels[" + str(sector) + "].png")
        plt.close()

    @staticmethod
    def compute_pixels_curves(tpf):
        masks = np.zeros(
            (tpf.shape[1] * tpf.shape[2], tpf.shape[1], tpf.shape[2]),
            dtype="bool",
        )
        for i in range(tpf.shape[1] * tpf.shape[2]):
            masks[i][np.unravel_index(i, (tpf.shape[1], tpf.shape[2]))] = True
        pixel_list = []
        for j in range(tpf.shape[1] * tpf.shape[2]):
            lc = tpf.to_lightcurve(aperture_mask=masks[j])
            lc = lc.remove_outliers(sigma_upper=3, sigma_lower=float('inf'))
            if len(lc.remove_nans().flux) == 0:
                pixel_list.append(None)
            else:
                pixel_list.append(lc)

    @staticmethod
    def vetting_field_of_view_single(fov_process_input):
        """
        Plots FOV for one sector data. To be called by a multiprocessing queue.
        :param fov_process_input: wrapper for the sector data
        """
        maglim = 6
        search_radius = 40
        try:
            tpf = fov_process_input.tpf_source.download(cutout_size=(CUTOUT_SIZE, CUTOUT_SIZE))
            row = tpf.row
            column = tpf.column
            plt.close()
            fig = plt.figure(figsize=(6.93, 5.5))
            gs = gridspec.GridSpec(1, 3, height_ratios=[1], width_ratios=[1, 0.05, 0.01])
            gs.update(left=0.05, right=0.95, bottom=0.12, top=0.95, wspace=0.01, hspace=0.03)
            ax1 = plt.subplot(gs[0, 0])
            # TPF plot
            mean_tpf = np.mean(tpf.flux.value, axis=0)
            nx, ny = np.shape(mean_tpf)
            norm = ImageNormalize(stretch=stretching.LogStretch())
            division = np.int(np.log10(np.nanmax(tpf.flux.value)))
            splot = plt.imshow(np.nanmean(tpf.flux, axis=0) / 10 ** division, norm=norm, cmap="viridis", \
                               extent=[column, column + ny, row, row + nx], origin='lower', zorder=0)
            aperture = fov_process_input.apertures[tpf.sector]
            aperture = aperture if isinstance(aperture, np.ndarray) else np.array(aperture)
            aperture_boolean = ApertureExtractor.from_pixels_to_boolean_mask(aperture, column, row, tpf.shape[2],
                                                                             tpf.shape[1])
            Watson.plot_tpf(tpf, tpf.sector, aperture_boolean, fov_process_input.save_dir)
            maskcolor = 'salmon'
            logging.info("    --> Using SHERLOCK aperture for sector %s...", tpf.sector)
            if aperture is not None:
                for pixels in aperture:
                    ax1.add_patch(patches.Rectangle((pixels[0], pixels[1]),
                                                    1, 1, color=maskcolor, fill=True, alpha=0.4))
                    ax1.add_patch(patches.Rectangle((pixels[0], pixels[1]),
                                                    1, 1, color=maskcolor, fill=False, alpha=1, lw=2))
            # Gaia sources
            gaia_id, mag = tpfplotter.get_gaia_data(fov_process_input.ra, fov_process_input.dec,
                                                    search_radius=search_radius)
            r, res = tpfplotter.add_gaia_figure_elements(tpf, magnitude_limit=mag + np.float(maglim), targ_mag=mag)
            x, y, gaiamags = r
            x, y, gaiamags = np.array(x) + 0.5, np.array(y) + 0.5, np.array(gaiamags)
            size = 128.0 / 2 ** ((gaiamags - mag))
            plt.scatter(x, y, s=size, c='red', alpha=0.6, edgecolor=None, zorder=10)
            # Gaia source for the target
            this = np.where(np.array(res['Source']) == int(gaia_id))[0]
            plt.scatter(x[this], y[this], marker='x', c='white', s=32, zorder=11)
            # Legend
            add = 0
            if np.int(maglim) % 2 != 0:
                add = 1
            maxmag = np.int(maglim) + add
            legend_mags = np.linspace(-2, maxmag, np.int((maxmag + 2) / 2 + 1))
            fake_sizes = mag + legend_mags  # np.array([mag-2,mag,mag+2,mag+5, mag+8])
            for f in fake_sizes:
                size = 128.0 / 2 ** ((f - mag))
                plt.scatter(0, 0, s=size, c='red', alpha=0.6, edgecolor=None, zorder=10,
                            label=r'$\Delta m=$ ' + str(np.int(f - mag)))
            ax1.legend(fancybox=True, framealpha=0.7)
            # Source labels
            dist = np.sqrt((x - x[this]) ** 2 + (y - y[this]) ** 2)
            dsort = np.argsort(dist)
            for d, elem in enumerate(dsort):
                if dist[elem] < 6:
                    plt.text(x[elem] + 0.1, y[elem] + 0.1, str(d + 1), color='white', zorder=100)
            # Orientation arrows
            tpfplotter.plot_orientation(tpf)
            # Labels and titles
            plt.xlim(column, column + ny)
            plt.ylim(row, row + nx)
            plt.xlabel('Pixel Column Number', fontsize=16)
            plt.ylabel('Pixel Row Number', fontsize=16)
            plt.title('Coordinates ' + fov_process_input.target_title + ' - Sector ' + str(tpf.sector),
                      fontsize=16)  # + ' - Camera '+str(tpf.camera))  #
            # Colorbar
            cbax = plt.subplot(gs[0, 1])  # Place it where it should be.
            pos1 = cbax.get_position()  # get the original position
            pos2 = [pos1.x0 - 0.05, pos1.y0, pos1.width, pos1.height]
            cbax.set_position(pos2)  # set a new position
            cb = Colorbar(ax=cbax, cmap="viridis", mappable=splot, orientation='vertical', ticklocation='right')
            plt.xticks(fontsize=14)
            exponent = r'$\times 10^' + str(division) + '$'
            cb.set_label(r'Flux ' + exponent + r' (e$^-$)', labelpad=10, fontsize=16)
            plt.savefig(fov_process_input.save_dir + '/fov_TPF_Gaia_' + fov_process_input.target_title + '_S' +
                        str(tpf.sector) + '.png')
            # Save Gaia sources info
            dist = np.sqrt((x - x[this]) ** 2 + (y - y[this]) ** 2)
            GaiaID = np.array(res['Source'])
            srt = np.argsort(dist)
            x, y, gaiamags, dist, GaiaID = x[srt], y[srt], gaiamags[srt], dist[srt], GaiaID[srt]
            IDs = np.arange(len(x)) + 1
            inside = np.zeros(len(x))
            for pixels in aperture:
                xtpf, ytpf = pixels[0], pixels[1]
                _inside = np.where((x > xtpf) & (x < xtpf + 1) &
                                   (y > ytpf) & (y < ytpf + 1))[0]
                inside[_inside] = 1
            data = Table([IDs, GaiaID, x, y, dist, dist * 21., gaiamags, inside.astype('int')],
                         names=['# ID', 'GaiaID', 'x', 'y', 'Dist_pix', 'Dist_arcsec', 'Gmag', 'InAper'])
            ascii.write(data, fov_process_input.save_dir + '/fov_Gaia_' + fov_process_input.target_title + '_S' +
                        str(tpf.sector) + '.dat', overwrite=True)
        except SystemExit:
            logging.exception("Field Of View generation tried to exit.")
        except Exception as e:
            logging.exception("Exception found when generating Field Of View plots")

    @staticmethod
    def vetting_field_of_view(indir, mission, tic, cadence, ra, dec, sectors, source, apertures,
                              cpus=multiprocessing.cpu_count() - 1):
        """
        Runs TPFPlotter to get field of view data.
        :param indir: the data source directory
        :param mission: the mission of the target
        :param tic: the target id
        :param cadence: the exposure time between measurements in seconds
        :param ra: the right ascension of the target
        :param dec: the declination of the target
        :param sectors: the sectors where the target was observed
        :param source: the source where the aperture was generated [tpf, tesscut, eleanor]
        :param apertures: a dict mapping sectors to boolean apertures
        :param cpus: cores to be used
        :return: the directory where resulting data is stored
        """
        try:
            sectors = [sectors] if isinstance(sectors, int) else sectors
            sectors_search = None if sectors is not None and len(sectors) == 0 else sectors
            logging.info("Preparing target pixel files for field of view plots")
            if mission != "TESS":
                return
            target_title = "TIC " + str(tic)
            #TODO use retrieval method depending on source parameter
            if cadence > 120:
                tpf_source = lightkurve.search_tesscut(target_title, sector=sectors_search)
                if tpf_source is None or len(tpf_source) == 0:
                    ra_str = str(ra)
                    dec_str = "+" + str(dec) if dec >= 0 else str(dec)
                    coords_str = ra_str + " " + dec_str
                    tpf_source = lightkurve.search_tesscut(coords_str, sector=sectors_search)
                    target_title = "RA={:.4f},DEC={:.4f}".format(ra, dec)
            else:
                tpf_source = lightkurve.search_targetpixelfile(target_title, sector=sectors_search,
                                                               author=lcbuilder.constants.SPOC_AUTHOR,
                                                               cadence=cadence)
            save_dir = indir
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            fov_process_inputs = []
            for i in range(0, len(tpf_source)):
                fov_process_inputs.append(FovProcessInput(save_dir, mission, tic, cadence, ra, dec, sectors, source,
                                                          apertures, tpf_source[i], target_title))
            with multiprocessing.Pool(processes=cpus) as pool:
                pool.map(Watson.vetting_field_of_view_single, fov_process_inputs)
            return save_dir
        except SystemExit:
            logging.exception("Field Of View generation tried to exit.")
        except Exception as e:
            logging.exception("Exception found when generating Field Of View plots")

class SingleTransitProcessInput:
    def __init__(self, data_dir, id, index, lc_file, lc_data_file, tpfs_dir, apertures,
                                         transit_times, depth, duration, period, rp_rstar, a_rstar, transits_mask):
        self.data_dir = data_dir
        self.id = id
        self.index = index
        self.lc_file = lc_file
        self.lc_data_file = lc_data_file
        self.tpfs_dir = tpfs_dir
        self.apertures = apertures
        self.transit_times = transit_times
        self.depth = depth
        self.duration = duration
        self.period = period
        self.rp_rstar = rp_rstar
        self.a_rstar = a_rstar
        self.transits_mask = transits_mask

class FovProcessInput:
    def __init__(self, save_dir, mission, tic, cadence, ra, dec, sectors, source, apertures, tpf_source, target_title):
        self.save_dir = save_dir
        self.mission = mission
        self.tic = tic
        self.cadence = cadence
        self.ra = ra
        self.dec = dec
        self.sectors = sectors
        self.source = source
        self.apertures = apertures
        self.tpf_source = tpf_source
        self.target_title = target_title
