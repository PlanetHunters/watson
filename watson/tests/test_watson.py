import multiprocessing
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch
from importlib import resources

import numpy
import pandas as pd
import pytest
import yaml
from lcbuilder import constants
from lcbuilder.star.starinfo import StarInfo

from watson.watson import TriceratopsThreadValidator, ValidatorInput, Watson


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

    @pytest.mark.integration
    def test_vetting_by_params(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199,
                                                    6.082, 0.25, [1, 2], 0.07571, cadence=[120],
                                                    author=[constants.SPOC_AUTHOR],
                                                    cpus=multiprocessing.cpu_count() // 2, clean=False)
            files_in_dir = os.listdir(vetting_dir)
            self.assertTrue({"transits_validation_report.pdf", "transits_validation_report_summary.pdf"}
                            .issubset(files_in_dir))
            self.assertEqual(len(files_in_dir), 33)
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    @pytest.mark.integration
    def test_vetting_by_files(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 8.32, 1327.51, 199, 6.082, 0.25,
                                                    [1, 2], 0.07571, author=[constants.SPOC_AUTHOR],
                                       a_rstar=20, cadence=[120], lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       star_file=object_dir + "params_star.csv",
                                       cpus=multiprocessing.cpu_count() // 2, clean=False, only_summary=True)
            files_in_dir = os.listdir(vetting_dir)
            self.assertTrue({"metrics.csv", "folded_cadences.png", "odd_even_folded_curves.png", "triceratops"}
                            .issubset(files_in_dir))
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    @pytest.mark.integration
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
                                                                      -69.226593, [1], "tpf", apertures,
                                                                      1)
                files_in_dir = os.listdir(fov_dir)
                self.assertGreaterEqual(len(files_in_dir), 1,
                    f"Expected at least 1 FOV file (TPF plot), got {len(files_in_dir)}. "
                    "The triceratops field plot may be missing if TRILEGAL.org is unreachable.")
        finally:
            if os.path.exists(fov_dir):
                shutil.rmtree(fov_dir, ignore_errors=False)

    @pytest.mark.integration
    def test_vetting_by_files_with_fov(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            with patch.object(Watson, "vetting_field_of_view") as fov_mock:
                def create_fov_outputs(indir, *args, **kwargs):
                    open(indir + "/fov_Flux_pixels[1].png", "wb").close()
                    open(indir + "/fov_sector_1.pdf", "wb").close()

                fov_mock.side_effect = create_fov_outputs
                Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199,
                                                        6.082, 0.25, [1, 2], 0.07571,
                                           a_rstar=20, cadence=[120], author=[constants.SPOC_AUTHOR],
                                                        lc_file=object_dir + "/lc.csv",
                                           lc_data_file=object_dir + "/lc_data.csv",
                                           tpfs_dir=object_dir + "/tpfs",
                                           apertures_file=object_dir + "/apertures.yaml",
                                           star_file=object_dir + "params_star.csv",
                                           cpus=multiprocessing.cpu_count() // 2, create_fov_plots=True,
                                           cadence_fov=120, ra=63.3739396231274, dec=-69.226822697583, clean=False)
                fov_mock.assert_called_once()
            files_in_dir = os.listdir(vetting_dir)
            self.assertTrue({"transits_validation_report.pdf", "transits_validation_report_summary.pdf"}
                            .issubset(files_in_dir))
            self.assertTrue(any(file.startswith("fov_Flux_pixels[") and file.endswith(".png")
                              for file in files_in_dir))
            self.assertTrue(any(file.startswith("fov_sector_") for file in files_in_dir))
        finally:
            if os.path.exists(vetting_dir):
                shutil.rmtree(vetting_dir, ignore_errors=False)

    @pytest.mark.integration
    def test_vetting_by_files_with_transits_list(self):
        object_dir = TestsWatson.get_path("TIC25155310_[1,_2]")
        vetting_dir = object_dir + "/vetting_0/"
        try:
            transits_list_df = pd.read_csv(object_dir + "/transits_stats.csv")
            transits_list_df = transits_list_df[transits_list_df["candidate"] == 0]
            Watson(object_dir, vetting_dir).vetting("TIC 25155310", 3.2899, 1327.51, 199,
                                                    6.082, 0.25, [1, 2], 0.07571,
                                                    author=[constants.SPOC_AUTHOR],
                                       a_rstar=20, cadence=[120], lc_file=object_dir + "/lc.csv",
                                       lc_data_file=object_dir + "/lc_data.csv",
                                       tpfs_dir=object_dir + "/tpfs",
                                       apertures_file=object_dir + "/apertures.yaml",
                                       star_file=object_dir + "params_star.csv",
                                       cpus=multiprocessing.cpu_count() // 2,
                                       transits_list=transits_list_df.to_dict("list"), ra=63.3739396231274,
                                        dec=-69.226822697583, clean=False)
            files_in_dir = os.listdir(vetting_dir)
            self.assertTrue({"transits_validation_report.pdf", "transits_validation_report_summary.pdf"}
                            .issubset(files_in_dir))
            self.assertTrue({"metrics.csv", "transit_depths.png", "triceratops"}.issubset(files_in_dir))
            self.assertTrue(any(file.startswith("single_transit_") and file.endswith(".png")
                               for file in files_in_dir))
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
        return str(resources.files("watson.tests") / path)


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
        """The calibrated score is applied directly to planet and FP scenarios."""
        df = self._make_scenarios(tp_prob=0.8, eb_prob=0.2)
        combined_fpp, combined_nfpp, updated_df = Watson._apply_bayesian_update(0.1, df)
        expected_fpp = (0.2 * 0.9) / (0.8 * 0.1 + 0.2 * 0.9)
        self.assertAlmostEqual(combined_fpp, expected_fpp, places=6)
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
        original_fpp = 0.19
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

        Setup: TP=0.95, EB=0.05, p_NN=0.9
        TP_unnorm = 0.95 * 0.9 = 0.855
        EB_unnorm = 0.05 * 0.1 = 0.005
        Z = 0.855 + 0.005 = 0.86
        combined_fpp = 0.005 / 0.86
        """
        df = self._make_scenarios(tp_prob=0.95, eb_prob=0.05)
        combined_fpp, _, _ = Watson._apply_bayesian_update(0.9, df)
        tp_u = 0.95 * 0.9
        eb_u = 0.05 * 0.1
        expected = eb_u / (tp_u + eb_u)
        self.assertAlmostEqual(combined_fpp, expected, places=10)


class TestsTriceratopsCapabilities(unittest.TestCase):
    def test_validator_input_preserves_extended_triceratops_options(self):
        validator_input = ValidatorInput(
            "results",
            object(),
            [],
            [],
            0.001,
            3.2,
            0.004,
            {},
            0,
            "contrast.csv",
            ignore_ebs=True,
            resolved_companion="yes",
            ignore_background_stars=True,
            bound_stars=[123],
            external_lc_files=["ground_r.csv"],
            filt_lcs=["r"],
        )

        self.assertTrue(validator_input.ignore_ebs)
        self.assertEqual("yes", validator_input.resolved_companion)
        self.assertTrue(validator_input.ignore_background_stars)
        self.assertEqual([123], validator_input.bound_stars)
        self.assertEqual(["ground_r.csv"], validator_input.external_lc_files)
        self.assertEqual(["r"], validator_input.filt_lcs)

    def test_validator_forwards_external_lightcurves_and_contrast_curve(self):
        class FakeTarget:
            def __init__(self):
                self.stars = pd.DataFrame({"ID": ["1"], "plx": [1.0]})
                self.probs = pd.DataFrame({
                    "ID": [1],
                    "scenario": ["TP"],
                    "prob": [1.0],
                })
                self.star_num = numpy.zeros(1)
                self.u1 = numpy.zeros(1)
                self.u2 = numpy.zeros(1)
                self.fluxratio_EB = numpy.zeros(1)
                self.fluxratio_comp = numpy.zeros(1)
                self.calc_probs_kwargs = None

            def calc_probs(self, **kwargs):
                self.calc_probs_kwargs = kwargs

            def plot_fits(self, **kwargs):
                return None

        with tempfile.TemporaryDirectory() as output_dir:
            fake_target = FakeTarget()
            validator_input = ValidatorInput(
                output_dir,
                fake_target,
                numpy.array([0.0, 1.0]),
                numpy.array([1.0, 1.0]),
                0.001,
                3.2,
                0.004,
                {},
                0,
                "contrast.csv",
                external_lc_files=["ground_r.csv"],
                filt_lcs=["r"],
            )
            TriceratopsThreadValidator.validate(validator_input)

        self.assertEqual("contrast.csv", fake_target.calc_probs_kwargs["contrast_curve_file"])
        self.assertEqual(["ground_r.csv"], fake_target.calc_probs_kwargs["external_lc_files"])
        self.assertEqual(["r"], fake_target.calc_probs_kwargs["filt_lcs"])


if __name__ == '__main__':
    unittest.main()
