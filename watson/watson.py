# from __future__ import print_function, absolute_import, division
import logging
import multiprocessing
import traceback

import PIL
import batman
import foldedleastsquares
import lcbuilder
import lightkurve
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import wotan
import yaml
from lcbuilder.constants import CUTOUT_SIZE
from lcbuilder.helper import LcbuilderHelper
from lcbuilder.lcbuilder_class import LcBuilder
from lcbuilder.objectinfo.MissionFfiIdObjectInfo import MissionFfiIdObjectInfo
from lcbuilder.objectinfo.MissionObjectInfo import MissionObjectInfo
from lcbuilder.photometry.aperture_extractor import ApertureExtractor
from lcbuilder.star.HabitabilityCalculator import HabitabilityCalculator
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
from math import ceil

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
            if not ffi:
                lc_build = lc_builder.build(MissionObjectInfo(id, sectors, cadence=cadence,
                                                              initial_transit_mask=transits_mask), self.data_dir)
            else:
                lc_build = lc_builder.build(MissionFfiIdObjectInfo(id, sectors, cadence=cadence,
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
        report = Report(self.data_dir, file_name, id, ra, dec, t0, period, duration, depth, transits_list, None,
                        v, j, h, k, with_tpfs)
        report.create_report()
        file_name = "transits_validation_report_summary.pdf"
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
        lc, lc_data, tpfs = Watson.initialize_lc_and_tpfs(id, lc_file, lc_data_file, tpfs_dir,
                                                          transits_mask=transits_mask)
        apertures = None
        if os.path.exists(apertures_file):
            apertures = yaml.load(open(apertures_file), yaml.SafeLoader)
            apertures = apertures["sectors"]
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
            num_of_transits = int(ceil(((last_time - t0) / period)))
            transit_lists = t0 + period * np.arange(0, num_of_transits)
            time_as_array = lc.time.value
            plot_range = duration / 3600 * 2
            transits_in_data = [time_as_array[(transit > time_as_array - plot_range) & (transit < time_as_array + plot_range)] for
                                transit in transit_lists]
            transit_t0s_list = transit_lists[[len(transits_in_data_set) > 0 for transits_in_data_set in transits_in_data]]
        self.plot_folded_curve(self.data_dir, id, lc, period, t0, duration, depth / 1000, rp_rstar, a_rstar)
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
            lc_data = Watson.normalize_lc_data(lc_data)
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
            for tpf_file in os.listdir(tpfs_dir):
                tpf = TessTargetPixelFile(tpfs_dir + "/" + tpf_file) if mission == lcbuilder.constants.MISSION_TESS else \
                    KeplerTargetPixelFile(tpfs_dir + "/" + tpf_file)
                tpfs.append(tpf)
        return lc, lc_data, tpfs

    @staticmethod
    def normalize_lc_data(lc_data):
        logging.info("Normalizing lc_data")
        time = lc_data["time"].to_numpy()
        dif = time[1:] - time[:-1]
        jumps = np.where(np.abs(dif) > 1)[0]
        jumps = np.append(jumps, len(lc_data))
        previous_jump_index = 0
        for jumpIndex in jumps:
            token = lc_data["centroids_x"][previous_jump_index:jumpIndex]
            lc_data["centroids_x"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            token = lc_data["centroids_y"][previous_jump_index:jumpIndex]
            lc_data["centroids_y"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            token = lc_data["motion_x"][previous_jump_index:jumpIndex]
            lc_data["motion_x"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            token = lc_data["motion_y"][previous_jump_index:jumpIndex]
            lc_data["motion_y"][previous_jump_index:jumpIndex] = token - np.nanmedian(token)
            previous_jump_index = jumpIndex
        return lc_data

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
        lc, lc_data, tpfs = Watson.initialize_lc_and_tpfs(single_transit_process_input.id,
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
        if lc_data is not None:
            zoom_lc_data = lc_data[(lc_data["time"] > transit_time - plot_range) & (lc_data["time"] < transit_time +
                                                                                plot_range)]
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
                tpf_short_framed.plot_pixels(aperture_mask=aperture_mask)
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
    def compute_phased_values_and_fill_plot(id, axs, lc, period, epoch, depth, duration, rp_rstar, a_rstar, range=5,
                                            bins=None):
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
            bin_stds, _, _ = stats.binned_statistic(folded_phase, folded_y, statistic='std', bins=bins)
            bin_nan_args = np.isnan(bin_stds)
            axs.errorbar(bin_centers[~bin_nan_args], bin_means[~bin_nan_args],
                         yerr=bin_stds[~bin_nan_args] / 2, xerr=bin_width / 2, marker='o', markersize=4,
                         color='darkorange', alpha=1, linestyle='none')
        else:
            bin_centers = folded_phase
            bin_means = folded_y
            bin_stds = folded_y_err * 2
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
        return axs, bin_centers, bin_means, bin_stds / 2

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
                tpf_source = lightkurve.search_targetpixelfile(target_title, sector=sectors_search, author="SPOC",
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
