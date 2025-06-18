import multiprocessing
import os
import shutil
import unittest

import numpy
import pandas as pd
import pkg_resources
import yaml
from lcbuilder import constants
from lcbuilder.star.starinfo import StarInfo

from watson.watson import Watson


class TestsWatson(unittest.TestCase):
    def test_iatson(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = TestsWatson.get_path("vetting_test")
        predictions_df, first_row_df, branches_results_df, values_results_df = Watson.run_iatson("TIC 25155310", 3.2899, 199, 1327.51, 6.082,
                          vetting_dir, object_dir + '/params_star.csv', object_dir + '/lc.csv', transits_mask=None, plot_inputs=False, batch_size=5)
        self.assertAlmostEqual(first_row_df['prediction_value_mean'].iloc[0], 0.2127, 3)
        self.assertAlmostEqual(first_row_df['prediction_value_std'].iloc[0], 0.1682, 3)
        self.assertAlmostEqual(first_row_df['prediction_value_cal_mean'].iloc[0], 0.17882, 3)
        self.assertAlmostEqual(first_row_df['prediction_value_cal_std'].iloc[0], 0.14567, 3)

    def test_vetting_by_params(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199,
                                                    6.082, [1, 2], 0.07571, cadence=120,
                                                    cpus=multiprocessing.cpu_count() // 2, clean=False)
            files_in_dir = os.listdir(vetting_dir)
            self.assertEqual(len(files_in_dir), 35)
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_vetting_by_files(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 8.32, 1327.51, 199, 6.082, [1, 2], 0.07571,
                                       a_rstar=20, cadence=120, lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       star_file=object_dir + "params_star.csv",
                                       cpus=multiprocessing.cpu_count() // 2, clean=False, only_summary=True)
            files_in_dir = os.listdir(vetting_dir)
            self.assertEqual(len(files_in_dir), 8)
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_fov_plots(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        fov_dir = object_dir + "/fov"
        os.mkdir(fov_dir)
        try:
            with open(object_dir + "/apertures.yaml") as apertures_file:
                apertures = yaml.load(apertures_file, yaml.SafeLoader)
                apertures = apertures["sectors"]
                Watson(object_dir, vetting_dir).vetting_field_of_view(fov_dir, "TESS", "25155310", 120, 63.374706,
                                                                      -69.226593, [1, 2], "tpf", apertures,
                                                                      1)
                files_in_dir = os.listdir(fov_dir)
                self.assertEqual(len(files_in_dir), 4)
        finally:
            if os.path.exists(fov_dir):
                shutil.rmtree(fov_dir, ignore_errors=False)

    def test_vetting_by_files_with_fov(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, [1, 2], 0.07571,
                                       a_rstar=20, cadence=120, lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       star_file=object_dir + "params_star.csv",
                                       cpus=multiprocessing.cpu_count() // 2, create_fov_plots=True,
                                       cadence_fov=120, ra=63.3739396231274, dec=-69.226822697583, clean=False)
            files_in_dir = os.listdir(vetting_dir)
            self.assertEqual(len(files_in_dir), 36)
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_vetting_by_files_with_transits_list(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            transits_list_df = pd.read_csv(object_dir + "/transits_stats.csv")
            transits_list_df = transits_list_df[transits_list_df["candidate"] == 0]
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, [1, 2], 0.07571,
                                       a_rstar=20, cadence=120, lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       star_file=object_dir + "params_star.csv",
                                       cpus=multiprocessing.cpu_count() // 2,
                                       transits_list=transits_list_df.to_dict("list"), ra=63.3739396231274,
                                       dec=-69.226822697583, clean=False)
            files_in_dir = os.listdir(vetting_dir)
            self.assertEquals(len(files_in_dir), 33)
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_create_report(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = TestsWatson.get_path("vetting_report_test/")
        transits_list_df = pd.read_csv(object_dir + "/transits_stats.csv")
        transits_list_df = transits_list_df[transits_list_df["candidate"] == 0]
        transits_list_df = transits_list_df[transits_list_df["depth"].notnull()]
        transits_list_df = transits_list_df[transits_list_df["t0"].notnull()]
        transits_list = numpy.array(transits_list_df.to_dict("list")["t0"])
        transits_list = transits_list[~numpy.isnan(transits_list)]
        Watson(vetting_dir, vetting_dir).report("TIC 25155310", 63.3739396231274, -69.226822697583, 1327.51, 3.2899, 199, 6.082,
                                   transits_list, [2, 5, 7], None, None, None, None)
        files_in_dir = os.listdir(vetting_dir)
        self.assertEqual(len(files_in_dir), 20)

    def test_compute_bootstrap_fap(self):
        vetting_dir = TestsWatson.get_path("vetting_report_test/")
        fap = (Watson(vetting_dir, vetting_dir)
         .compute_bootstrap_fap(numpy.linspace(0, 100, 100),
                                numpy.ones(100), 0.05, 0.05,
                                StarInfo(radius=1, mass=1)))
        self.assertAlmostEqual(fap, 0.00, places=2)

    @staticmethod
    def get_path(path):
        """
        Gets right path for tests environment
        :param path:
        :return: the real path of the test resource
        """
        return pkg_resources.resource_filename(__name__, path)


if __name__ == '__main__':
    unittest.main()
