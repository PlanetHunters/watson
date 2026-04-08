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
        self.assertAlmostEqual(first_row_df['prediction_value_mean'].iloc[0], 0.9015, 3)
        self.assertAlmostEqual(first_row_df['prediction_value_std'].iloc[0], 0.13111425171600377, 3)
        self.assertAlmostEqual(first_row_df['prediction_value_cal_mean'].iloc[0], 0.996926236152649, 3)
        self.assertAlmostEqual(first_row_df['prediction_value_cal_std'].iloc[0], 0.00540493406486552, 3)

    def test_vetting_by_params(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199,
                                                    6.082, 0.25, [1, 2], 0.07571, cadence=120,
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
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 8.32, 1327.51, 199, 6.082, 0.25, [1, 2], 0.07571,
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
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, 0.25, [1, 2], 0.07571,
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
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199, 6.082, 0.25, [1, 2], 0.07571,
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


class TestsBayesianFpp(unittest.TestCase):
    """Unit tests for Watson._apply_bayesian_update and Watson.compute_bayesian_fpp."""

    def _make_scenarios(self, tp_prob, eb_prob, stp_prob=0.0, neb_prob=0.0):
        """Build a minimal scenarios DataFrame with known raw probabilities."""
        rows = [
            {'scenario': 'TP', 'prob': tp_prob},
            {'scenario': 'EB', 'prob': eb_prob},
        ]
        if stp_prob > 0:
            rows.append({'scenario': 'STP', 'prob': stp_prob})
        if neb_prob > 0:
            rows.append({'scenario': 'NEB', 'prob': neb_prob})
        return pd.DataFrame(rows)

    def test_neutral_score(self):
        """p_NN = pi_train should leave FPP unchanged (BF = 1, no information added)."""
        df = self._make_scenarios(tp_prob=0.8, eb_prob=0.2)
        original_fpp = 0.2  # EB fraction
        pi_train = 0.1
        # neutral point: p_NN == pi_train → planet_likelihood == fp_likelihood → no update
        combined_fpp, combined_nfpp, updated_df = Watson._apply_bayesian_update(pi_train, df, pi_train)
        self.assertAlmostEqual(combined_fpp, original_fpp, places=6)
        self.assertAlmostEqual(combined_nfpp, 0.0, places=6)

    def test_high_score_reduces_fpp(self):
        """High p_NN (planet-like signal) should lower FPP."""
        df = self._make_scenarios(tp_prob=0.8, eb_prob=0.2)
        original_fpp = 0.2
        combined_fpp, _, _ = Watson._apply_bayesian_update(0.9, df)
        self.assertLess(combined_fpp, original_fpp)

    def test_low_score_increases_fpp(self):
        """Low p_NN (EB-like signal) should raise FPP."""
        df = self._make_scenarios(tp_prob=0.8, eb_prob=0.2)
        original_fpp = 0.2
        combined_fpp, _, _ = Watson._apply_bayesian_update(0.1, df)
        self.assertGreater(combined_fpp, original_fpp)

    def test_off_target_also_reduced_with_high_score(self):
        """Off-target scenarios (STP, NEB) should also decrease when p_NN is high."""
        df = self._make_scenarios(tp_prob=0.6, eb_prob=0.2, stp_prob=0.1, neb_prob=0.1)
        _, _, updated_df = Watson._apply_bayesian_update(0.9, df)
        # STP prior = 0.1/1.0 = 0.1 of total; after update it should be smaller
        stp_after = updated_df.loc[updated_df['scenario'] == 'STP', 'prob'].values[0]
        neb_after = updated_df.loc[updated_df['scenario'] == 'NEB', 'prob'].values[0]
        # TP prior = 0.6 → after update TP fraction should be larger
        tp_after = updated_df.loc[updated_df['scenario'] == 'TP', 'prob'].values[0]
        self.assertGreater(tp_after, 0.6)
        self.assertLess(stp_after, 0.1)
        self.assertLess(neb_after, 0.1)

    def test_normalization(self):
        """Output scenario probabilities must always sum to 1."""
        df = self._make_scenarios(tp_prob=0.6, eb_prob=0.2, stp_prob=0.1, neb_prob=0.1)
        for nn_score in [0.05, 0.5, 0.95]:
            _, _, updated_df = Watson._apply_bayesian_update(nn_score, df)
            self.assertAlmostEqual(updated_df['prob'].sum(), 1.0, places=10)

    def test_missing_files_skips(self):
        """compute_bayesian_fpp() should return None silently if inputs are missing."""
        result = Watson.compute_bayesian_fpp('/nonexistent/path')
        self.assertIsNone(result)

    def test_known_values(self):
        """Verify exact combined FPP against manually computed value.

        Setup: TP=0.95, EB=0.05, p_NN=0.9, pi_train=0.1
        planet_likelihood = 0.9 / 0.1 = 9.0
        fp_likelihood     = 0.1 / 0.9 ≈ 0.1111
        TP_unnorm = 0.95 * 9.0   = 8.55
        EB_unnorm = 0.05 * 0.1111 ≈ 0.005556
        Z = 8.55 + 0.005556 = 8.555556
        combined_fpp = EB_unnorm / Z = 0.005556 / 8.555556
        """
        df = self._make_scenarios(tp_prob=0.95, eb_prob=0.05)
        combined_fpp, _, _ = Watson._apply_bayesian_update(0.9, df, training_planet_fraction=0.1)
        tp_u = 0.95 * (0.9 / 0.1)
        eb_u = 0.05 * (0.1 / 0.9)
        expected = eb_u / (tp_u + eb_u)
        self.assertAlmostEqual(combined_fpp, expected, places=10)


if __name__ == '__main__':
    unittest.main()
