import os
import re

import lcbuilder.eleanor
import sys
sys.modules['eleanor'] = sys.modules['lcbuilder.eleanor']
import eleanor
from lcbuilder.eleanor.targetdata import TargetData
from lcbuilder.eleanor_manager import EleanorManager

import foldedleastsquares
import lcbuilder
import lightkurve
import pandas as pd
import logging
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs
from astropy import units as u
from lcbuilder.constants import LIGHTKURVE_CACHE_DIR, ELEANOR_CACHE_DIR, CUTOUT_SIZE
from lcbuilder.objectinfo.preparer.MissionLightcurveBuilder import MissionLightcurveBuilder
from scipy.integrate import dblquad


def create_star_csv(create_star_input):
    tries = 0
    object_id = None
    mission_id = None
    while tries < 3:
        try:
            mission, mission_prefix, id_int = MissionLightcurveBuilder().parse_object_id(create_star_input.id)
            star_info = create_star_input.star_catalog.catalog_info(id_int)
            star_df = pd.DataFrame(
                columns=['id', 'ra', 'dec', 'ld_a', 'ld_b', 'Teff', 'lum', 'logg', 'radius', 'mass', 'v', 'j', 'h',
                         'k',
                         'dist_arcsec'])
            star_df = star_df.append({'id': object_id, 'ra': star_info[11], 'dec': star_info[12],
                                      'ld_a': star_info[0][0],
                                      'ld_b': star_info[0][1], 'Teff': star_info[1],
                                      'lum': star_info[2], 'logg': star_info[3], 'radius': star_info[5],
                                      'mass': star_info[8], 'v': star_info[13], 'j': star_info[15],
                                      'h': star_info[17],
                                      'k': star_info[19], 'dist_arcsec': 0}, ignore_index=True)
            break
        except Exception as e:
            logging.exception("Failed object %s try", object_id)
            tries = tries + 1
    if tries >= 3:
        print("Failed downloading object id " + str(object_id))
    assert tries < 3
    try:  # Acquiring neighbours parameters
        ra = star_df['ra'].iloc[0]
        dec = star_df['dec'].iloc[0]
        if ra is not None and dec is not None:
            ticid = Catalogs.query_region(
                SkyCoord(ra, dec, unit="deg"),
                radius=create_star_input.search_radius * create_star_input.pixel_size,
                catalog="TIC"
            )[0]["ID"]
        df = Catalogs.query_object(
            "TIC" + str(ticid),
            radius=create_star_input.search_radius * create_star_input.pixel_size,
            catalog="TIC"
        )
        stars = df.to_pandas()
        sep = [0]
        pa = [0]
        c_target = SkyCoord(
            stars["ra"].values[0],
            stars["dec"].values[0],
            unit="deg"
        )
        for i in range(1, len(stars)):
            c_star = SkyCoord(
                stars["ra"].values[i],
                stars["dec"].values[i],
                unit="deg"
            )
            sep.append(
                np.round(
                    c_target.separation(c_star).to(u.arcsec).value,
                    3
                )
            )
            pa.append(
                np.round(
                    c_target.position_angle(c_star).to(u.deg).value,
                    3
                )
            )
        stars["dist_arcsec"] = sep
        stars["PA (E of N)"] = pa  # TODO should we use this?
        for index, star in stars.iterrows():
            star_df = star_df.append({'id': star['ID'], 'ra': star['ra'], 'dec': star['dec'], 'ld_a': 0,
                                      'ld_b': 0, 'Teff': star['Teff'],
                                      'lum': star['lum'], 'logg': star['logg'], 'radius': star['rad'],
                                      'mass': star['mass'], 'v': star['Vmag'], 'j': star['Jmag'],
                                      'h': star['Hmag'],
                                      'k': star['Kmag'], 'dist_arcsec': star['dist_arcsec']}, ignore_index=True)
        star_df = star_df.sort_values(["dist_arcsec"], ascending=True)
    except Exception as e:
        logging.exception('Something failed when retrieving neighbours info for object id %s', str(object_id))
    file = create_star_input.output_dir + '/' + str(create_star_input.id) + '_' + str(object_id) + '_stars.csv'
    star_df.to_csv(file)
    logging.info('Processed object id ' + str(object_id))
    return file


def get_neighbour_lc(neighbour_input):
    lc_df = pd.DataFrame(columns=["folded_time", "flux"])
    lc = None
    if neighbour_input.mission == "TESS":
        star_lcs = lightkurve.search_lightcurve(
            "TIC " + str(neighbour_input.id),
            mission=neighbour_input.mission,
            sector=None,
            author=neighbour_input.author,
            cadence="long"
        ).download_all(download_dir=os.path.expanduser('~') + '/' + LIGHTKURVE_CACHE_DIR)
        if star_lcs is None:
            star = eleanor.multi_sectors(tic=round(neighbour_input.id), sectors='all',
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
            #return lc_df
        else:
            matching_objects = []
            for i in range(0, len(star_lcs.data)):
                if star_lcs.data[i].label == "TIC " + str(neighbour_input.id):
                    if lc is None:
                        lc = star_lcs.data[i].normalize()
                    else:
                        lc = lc.append(star_lcs.data[i].normalize())
                else:
                    matching_objects.append(star_lcs.data[i].label)
            if lc is None:
                return lc_df
            lc = lc.remove_nans()
            if neighbour_input.mission == lcbuilder.constants.MISSION_K2:
                lc = lc.to_corrector("sff").correct(windows=20)
    else:
        star_lcs = lightkurve.search_lightcurve(
            "TIC " + str(neighbour_input.id),
            mission=neighbour_input.mission,
            sector=None,
            author=neighbour_input.author,
            cadence="long"
        ).download_all(download_dir=os.path.expanduser('~') + '/' + LIGHTKURVE_CACHE_DIR)
        if star_lcs is None:
            return lc_df
        matching_objects = []
        for i in range(0, len(star_lcs.data)):
            if star_lcs.data[i].label == "TIC " + str(neighbour_input.id):
                if lc is None:
                    lc = star_lcs.data[i].normalize()
                else:
                    lc = lc.append(star_lcs.data[i].normalize())
            else:
                matching_objects.append(star_lcs.data[i].label)
        if lc is None:
            return lc_df
        lc = lc.remove_nans()
        if neighbour_input.mission == lcbuilder.constants.MISSION_K2:
            lc = lc.to_corrector("sff").correct(windows=20)
    fold_times = foldedleastsquares.fold(lc.time.value, neighbour_input.period,
                                         neighbour_input.epoch + neighbour_input.period / 2)
    lc_df['folded_time'] = fold_times
    lc_df['flux'] = lc.flux.value
    lc_df = lc_df.sort_values(by=['folded_time'])
    lc_df = lc_df[(lc_df['folded_time'] > 0.5 - 3 * neighbour_input.duration) &
                  (lc_df['folded_time'] < 0.5 + 3 * neighbour_input.duration)]
    return lc_df

class NeighbourInput:

    def __init__(self, index, mission, author, id, period, epoch, duration) -> None:
        self.index = index
        self.mission = mission
        self.author = author
        self.id = id
        self.period = period
        self.epoch = epoch
        self.duration = duration


# def calc_depths(tdepth: float, stars_file: str, all_ap_pixels=None):
#     """
#     Calculates the transit depth each source near the target would
#     have if it were the source of the transit.
#     This is done by modeling the PSF of each source as a circular
#     Gaussian with a standard deviation of 0.75 pixels.
#     Args:
#         tdepth (float): Reported transit depth [ppm].
#         all_ap_pixels (list of numpy arrays): Apertures used to
#                                               extract light curve.
#     """
#     assert all_ap_pixels is not None
#     assert stars_file is not None
#     stars = pd.read_csv(stars_file)
#     # for each aperture, calculate contribution due to each star
#     rel_flux = np.zeros([len(stars)])
#     flux_ratio = np.zeros([len(stars)])
#     largest_aperture = []
#     aperture_sector = None
#     for sector, ap_pixels in all_ap_pixels:
#         if len(np.argwhere(ap_pixels is True).flatten()) > len(np.argwhere(largest_aperture is True).flatten()):
#             largest_aperture = ap_pixels
#             aperture_sector = sector
#     for i in range(len(stars)):
#         star_lcs = lightkurve.search_lightcurve(
#             stars.loc[i, "id"],
#             mission="TESS",
#             sectors="all",
#             cadence="long"
#         ).download_all(download_dir=LIGHTKURVE_CACHE_DIR)
#         star_tpfs = lightkurve.search_targetpixelfile(
#             stars.loc[i, "id"],
#             mission="TESS",
#             sectors="all",
#             cadence="long"
#         ).download_all(download_dir=LIGHTKURVE_CACHE_DIR)
#         for star_tpf in star_tpfs:
#             if star_tpf.sector == aperture_sector:
#
#         # location of star in pixel space for aperture k
#         mu_x = self.pix_coords[k][i, 0]
#         mu_y = self.pix_coords[k][i, 1]
#         # star's flux normalized to brightest star
#         A = 10 ** ((np.min(stars.v.values) - stars.v.values[i]) / 2.5)
#         # integrate PSF in each pixel
#         this_flux = 0
#         for j in range(len(largest_aperture)):
#             this_pixel = largest_aperture[j]
#             this_flux += dblquad(
#                 Gauss2D,
#                 this_pixel[1] - 0.5,
#                 this_pixel[1] + 0.5,
#                 this_pixel[0] - 0.5,
#                 this_pixel[0] + 0.5,
#                 args=(mu_x, mu_y, 0.75, A))[0]
#         rel_flux[i] = this_flux
#         # calculate flux ratios for this aperture
#         flux_ratio[:] = (rel_flux[:] / np.sum(rel_flux))
#
#     # take average of flux ratios across all apertures and append
#     # to stars dataframe
#     stars["fluxratio"] = flux_ratio
#     # calculate transit depth of each star given input transit depth
#     tdepths = np.zeros(len(stars))
#     for i in range(len(flux_ratio)):
#         if flux_ratio[i] != 0:
#             tdepths[i] = 1 - (flux_ratio[i] - tdepth) / flux_ratio[i]
#     tdepths[tdepths > 1] = 0
#     stars["tdepth"] = tdepths
#
#     # check target and possible NFPs for missing properties
#     filtered_stars = stars[stars["tdepth"] > 0]
#     for i, ID in enumerate(filtered_stars["ID"].values):
#         M_s = filtered_stars["mass"].values[i]
#         R_s = filtered_stars["rad"].values[i]
#         Teff = filtered_stars["Teff"].values[i]
#         Tmag = filtered_stars["Tmag"].values[i]
#         plx = filtered_stars["plx"].values[i]
#         if i == 0:
#             if (np.isnan(M_s) or np.isnan(R_s)
#                     or np.isnan(Teff) or np.isnan(plx)):
#                 print(
#                     "WARNING: " + str(ID)
#                     + " is missing stellar properties required "
#                     + "for validation."
#                     + "Please ensure a stellar "
#                     + "mass (in M_Sun), radius (in R_Sun), "
#                     + "Teff (in K), and plx (in mas) "
#                     + "are provided in the .stars dataframe."
#                 )
#         else:
#             if (np.isnan(M_s) or np.isnan(R_s)
#                     or np.isnan(Teff)):
#                 print(
#                     "WARNING: " + str(ID)
#                     + " is missing stellar properties. "
#                     + "If a mass (in M_Sun), "
#                     + "radius (in R_Sun), and/or Teff (in K) "
#                     + "are not added to the .stars dataframe, "
#                     + "Solar values will be assumed."
#                 )
#     return


def Gauss2D(x, y, mu_x, mu_y, sigma, A):
    """
    Calculates a circular Gaussian at specified grid points.
    Args:
        x, y (1D numpy arrays): Grid that you would like to calculate
                                Gaussian over.
        mu_x, mu_y (floats): Locations of star / Gaussian peak.
        sigma (float): Standard deviation of Gaussian
        A (float): Area under Gaussian.
    Returns:
    """
    xgrid, ygrid = np.meshgrid(x, y)
    exponent = ((xgrid - mu_x) ** 2 + (ygrid - mu_y) ** 2) / (2 * sigma ** 2)
    GaussExp = np.exp(-exponent)
    return A / (2 * np.pi * sigma ** 2) * GaussExp


class CreateStarCsvInput:
    def __init__(self, lc_file, mission, id, pixel_size, search_radius, lcs_regex, star_catalog, output_dir) -> None:
        super().__init__()
        self.lc_file = lc_file
        self.mission = mission
        self.id = id
        self.search_radius = search_radius
        self.pixel_size = pixel_size * u.arcsec
        self.lcs_regex = lcs_regex
        self.star_catalog = star_catalog
        self.output_dir = output_dir
