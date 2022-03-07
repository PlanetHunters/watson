import multiprocessing
import os
import shutil
import unittest
import pandas as pd
import pkg_resources
import yaml

from patoso.patoso import Patoso


class TestsPatoso(unittest.TestCase):
    def test_vetting_by_params(self):
        object_dir = self.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0"
        try:
            Patoso(object_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, False, [1, 2], 0.07571,
                                       cadence=120, cpus=multiprocessing.cpu_count() // 2)
            files_in_dir = os.listdir(vetting_dir)
            assert len(files_in_dir) == 17
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_vetting_by_files(self):
        object_dir = self.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0"
        try:
            Patoso(object_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, False, [1, 2], 0.07571,
                                       a_rstar=20, cadence=120, lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       cpus=multiprocessing.cpu_count() // 2)
            files_in_dir = os.listdir(vetting_dir)
            assert len(files_in_dir) == 17
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_fov_plots(self):
        object_dir = self.get_path("TIC25155310_[1,_2]")
        fov_dir = object_dir + "/fov"
        os.mkdir(fov_dir)
        try:
            apertures = yaml.load(open(object_dir + "/apertures.yaml"), yaml.SafeLoader)
            apertures = apertures["sectors"]
            Patoso(object_dir).vetting_field_of_view(fov_dir, "TESS", "25155310", 120, 63.374706, -69.226593, [1, 2],
                                                     "tpf", apertures, multiprocessing.cpu_count() // 2)
            files_in_dir = os.listdir(fov_dir)
            assert len(files_in_dir) == 6
        finally:
            if os.path.exists(fov_dir):
                shutil.rmtree(fov_dir, ignore_errors=False)

    def test_vetting_by_params_with_fov(self):
        object_dir = self.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0"
        try:
            Patoso(object_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, False, [1, 2], 0.07571,
                                       cadence=120, cpus=multiprocessing.cpu_count() // 2, create_fov_plots=True,
                                       cadence_fov=1800, ra_fov=63.3739396231274, dec_fov=-69.226822697583)
            files_in_dir = os.listdir(vetting_dir)
            assert len(files_in_dir) == 23
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_vetting_by_files_with_fov(self):
        object_dir = self.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0"
        try:
            Patoso(object_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, False, [1, 2], 0.07571,
                                       a_rstar=20, cadence=120, lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       cpus=multiprocessing.cpu_count() // 2, create_fov_plots=True,
                                       cadence_fov=120, ra_fov=63.3739396231274, dec_fov=-69.226822697583)
            files_in_dir = os.listdir(vetting_dir)
            assert len(files_in_dir) == 23
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def test_vetting_by_files_with_transits_list(self):
        object_dir = self.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0"
        try:
            transits_list_df = pd.read_csv(object_dir + "/transits_stats.csv")
            transits_list_df = transits_list_df[transits_list_df["candidate"] == 0]
            Patoso(object_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, False, [1, 2], 0.07571,
                                       a_rstar=20, cadence=120, lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       cpus=multiprocessing.cpu_count() // 2,
                                       transits_list=transits_list_df.to_dict("list"))
            files_in_dir = os.listdir(vetting_dir)
            assert len(files_in_dir) == 18
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    def get_path(self, path):
        return pkg_resources.resource_filename(__name__, path)


if __name__ == '__main__':
    unittest.main()
