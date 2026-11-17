"""The compact landscape vetting report (synthetic data: no network, no TPFs)."""
import logging
import os
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from watson import landscape_report as lr
from watson.tests.vetting_fixtures import (make_landscape_dir, make_vetting_dir, pdf_text, PERIOD, T0, DURATION_MIN,
                                           DEPTH_PPT)


def inputs_for(directory, lc_file=None, transits=None, **kw):
    return lr.ReportInputs(directory, "TIC 1", 63.37, -69.22, T0, PERIOD, DURATION_MIN, DEPTH_PPT, 12.8, 9.4, 8.8, 8.6,
                           lc_file=lc_file, transits=transits, **kw)


class NumericsTest(unittest.TestCase):
    def test_binned_ignores_sparse_bins_and_nan(self):
        x = np.array([0.1, 0.2, 0.3, 5.0])
        y = np.array([1.0, 3.0, np.nan, 7.0])
        cx, cy, _ = lr.binned(x, y, 0, 10, 2, minimum=2)
        self.assertEqual(len(cx), 1)
        self.assertAlmostEqual(cy[0], 2.0)

    def test_box_model_depth_and_edges(self):
        m = lr.box_model([0, 0.2, 5], 0.001, 1.0)
        self.assertAlmostEqual(m[0], 0.999)
        self.assertAlmostEqual(m[2], 1.0)
        self.assertTrue(0.999 < lr.box_model([0.45], 0.001, 1.0)[0] < 1.0)

    def test_phase_hours_wraps_to_mid_transit_and_secondary(self):
        self.assertAlmostEqual(float(lr.phase_hours(T0 + 5 * PERIOD + 1 / 24, T0, PERIOD)), 1.0, 6)
        self.assertAlmostEqual(float(lr.phase_hours(T0 + PERIOD / 2, T0, PERIOD, shift=0.5)), 0.0, 6)


class MetricsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_classifier_thresholds_match_the_old_report(self):
        make_vetting_dir(self.tmp.name)
        m = lr.classifier_metrics(self.tmp.name).set_index("metric")
        self.assertEqual(m.loc["WATSON-NET", "passed"], True)       # 0.994 + 0.011 >= 0.646
        self.assertEqual(m.loc["WATSON-NET err", "passed"], True)   # 0.011 < 0.015
        self.assertTrue(np.isnan(m.loc["Triceratops_FPP", "passed"]))  # 0.015 < 0.072 <= 0.5
        self.assertEqual(m.loc["Triceratops_NFPP", "passed"], True)
        self.assertEqual(m.loc["Bayesian_FPP", "passed"], True)

    def test_a_low_network_score_fails(self):
        make_vetting_dir(self.tmp.name)
        pd.DataFrame({"prediction_value_cal_mean": [0.01], "prediction_value_cal_std": [0.02]}) \
            .to_csv(f"{self.tmp.name}/iatson_averages.csv", index=False)
        m = lr.classifier_metrics(self.tmp.name).set_index("metric")
        self.assertEqual(m.loc["WATSON-NET", "passed"], False)
        self.assertTrue(np.isnan(m.loc["WATSON-NET err", "passed"]))  # 0.015 <= 0.02 < 0.1

    def test_load_metrics_and_status(self):
        make_vetting_dir(self.tmp.name)
        m = lr.load_metrics(self.tmp.name)
        self.assertEqual(lr.status_of(m, "snr"), "pass")
        self.assertEqual(lr.status_of(m, "secondary_snr"), "fail")
        self.assertEqual(lr.status_of(m, "odd_even_depth_sigma"), "none")
        self.assertEqual(lr.status_of(m, "does_not_exist"), "none")

    def test_verdict_names_the_failures_and_never_claims_a_planet(self):
        make_vetting_dir(self.tmp.name)
        headline, _, detail = lr.verdict(lr.load_metrics(self.tmp.name))
        self.assertIn("failed", headline)
        self.assertIn("secondary_snr", detail)
        ok = pd.DataFrame({"score": [1.0], "passed": [1]}, index=pd.Index(["snr"], name="metric"))
        self.assertEqual(lr.verdict(ok)[0], "No test failed")
        self.assertEqual(lr.verdict(ok.iloc[0:0])[0], "No verdict available")

    def test_empty_directory_gives_empty_metrics(self):
        self.assertEqual(len(lr.load_metrics(self.tmp.name)), 0)


class FlaggedTransitsTest(unittest.TestCase):
    def test_extremes_closest_and_outliers_are_flagged(self):
        transits = {"t0": list(range(10)), "depth": [6.0, 6.1, 5.9, 8.5, 6.0, 6.05, 5.95, 6.0, 3.0, 6.02],
                    "depth_err": [0.3] * 10}
        frame, z, flagged = lr.transit_depths_summary(transits, 6.0)
        self.assertEqual(len(frame), 10)
        self.assertTrue(flagged[3].startswith("deepest"))
        self.assertIn("outlier", flagged[3])
        self.assertIn(3, flagged)
        self.assertIn(8, flagged)
        self.assertTrue(any("outlier" in r for r in flagged.values()))
        self.assertLessEqual(len(flagged), 10)

    def test_nan_and_zero_error_transits_are_dropped(self):
        transits = {"t0": [1, 2, 3, 4], "depth": [6.0, np.nan, 6.1, 6.2], "depth_err": [0.3, 0.3, 0.0, 0.3]}
        frame, _, _ = lr.transit_depths_summary(transits, 6.0)
        self.assertEqual(frame.t0.tolist(), [1, 4])

    def test_empty_input(self):
        frame, z, flagged = lr.transit_depths_summary({"t0": [], "depth": [], "depth_err": []}, 6.0)
        self.assertEqual((len(frame), len(z), flagged), (0, 0, {}))


class WriteReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.lc_file, self.transits = make_landscape_dir(self.tmp.name)
        self.path = f"{self.tmp.name}/report.pdf"

    def build(self, **kw):
        return lr.write_landscape_report(inputs_for(self.tmp.name, self.lc_file, self.transits, **kw), self.path)

    def test_is_landscape_short_and_has_no_partial_file(self):
        pages = self.build()
        self.assertFalse(os.path.exists(self.path + ".partial"))
        n, text, (w, h) = pdf_text(self.path)
        self.assertEqual(n, pages)
        self.assertGreater(w, h)
        self.assertLessEqual(pages, 8)
        for title in ("Verdict and signal", "Classifiers", "Is the signal on the target?", "Single transits",
                      "Detailed sheets"):
            self.assertIn(title, text)

    def test_flagged_transit_gets_a_sheet_and_unflagged_does_not(self):
        self.build()
        _, text, _ = pdf_text(self.path)
        self.assertIn("Transit #3 ", text)         # the injected outlier
        self.assertIn("depth outlier", text)

    def test_verdict_page_lists_the_failing_metric(self):
        self.build()
        _, text, _ = pdf_text(self.path)
        self.assertIn("failed", text)
        self.assertIn("secondary_snr", text)

    def test_sector_pixel_images_are_shown_on_the_target_page(self):
        from PIL import Image
        for sector in (15, 22):
            Image.new("RGB", (280, 440), (120, 160, 200)).save(f"{self.tmp.name}/folded_tpf_{sector}_0.png")
        self.build()
        _, text, _ = pdf_text(self.path)
        self.assertIn("sector 15", text)
        self.assertIn("sector 22", text)

    def test_without_transits_the_transit_pages_are_skipped(self):
        lr.write_landscape_report(inputs_for(self.tmp.name, self.lc_file, None), self.path)
        n, text, _ = pdf_text(self.path)
        self.assertEqual(n, 3)
        self.assertNotIn("Single transits", text)

    def test_without_any_input_file_it_still_writes_a_report(self):
        with tempfile.TemporaryDirectory() as empty:
            path = f"{empty}/r.pdf"
            lr.write_landscape_report(inputs_for(empty), path)
            n, text, _ = pdf_text(path)
            self.assertEqual(n, 3)
            self.assertIn("No verdict available", text)

    def test_many_flagged_transits_are_capped(self):
        transits = {"t0": [T0 + i * PERIOD for i in range(60)],
                    "depth": [DEPTH_PPT + (3 if i % 2 else -3) for i in range(60)], "depth_err": [0.2] * 60}
        lr.write_landscape_report(inputs_for(self.tmp.name, self.lc_file, transits), self.path)
        n, text, _ = pdf_text(self.path)
        self.assertLessEqual(n, 4 + lr.MAX_SHEETS // 2)
        self.assertIn("more flagged transits without a sheet", text)

    def test_a_failing_page_becomes_a_stub_and_the_rest_survives(self):
        with mock.patch.object(lr, "page_classifiers", side_effect=RuntimeError("boom")), \
                self.assertLogs(level=logging.ERROR):
            self.build()
        n, text, _ = pdf_text(self.path)
        self.assertIn("Not available: RuntimeError: boom", text)
        self.assertIn("Verdict and signal", text)
        self.assertIn("Single transits", text)

    def test_pages_are_drawn_with_the_report_style_whatever_the_global_style(self):
        seen = []
        original = lr.PdfPages.savefig

        def spy(pdf, fig, *args, **kwargs):
            params = lr.matplotlib.rcParams
            seen.append((params["font.size"], params["axes.labelsize"], params["xtick.labelsize"],
                         params["ytick.labelsize"], params["legend.fontsize"]))
            return original(pdf, fig, *args, **kwargs)

        # what seaborn does: absolute sizes that ignore font.size
        hostile = {"font.size": 18, "axes.labelsize": 20, "xtick.labelsize": 16, "ytick.labelsize": 16,
                   "legend.fontsize": 15}
        with lr.matplotlib.rc_context(hostile), mock.patch.object(lr.PdfPages, "savefig", spy):
            self.build()
        self.assertTrue(seen)
        self.assertEqual(set(seen), {(lr.RC["font.size"], "medium", "medium", "medium", "medium")})

    def test_the_report_style_does_not_depend_on_the_current_rcparams(self):
        with lr.matplotlib.rc_context({"axes.labelsize": 30, "font.family": "serif"}):
            style = lr.report_style()
        self.assertEqual(style["axes.labelsize"], "medium")
        self.assertEqual(style["font.size"], lr.RC["font.size"])
        self.assertNotIn("backend", style)

    def test_the_report_is_far_smaller_than_the_old_one(self):
        self.build()
        self.assertLess(os.path.getsize(self.path), 5_000_000)


class WatsonWiringTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.lc_file, self.transits = make_landscape_dir(self.tmp.name)
        from watson.watson import Watson
        self.watson = Watson(self.tmp.name, self.tmp.name)
        self.call = lambda: self.watson._Watson__landscape_report(
            "TIC 1", 63.37, -69.22, T0, PERIOD, DURATION_MIN, DEPTH_PPT, self.transits, 12.8, 9.4, 8.8, 8.6,
            self.lc_file, self.tmp.name + "/none.csv", self.tmp.name + "/none", self.tmp.name + "/none.yaml",
            self.tmp.name + "/none.csv", None)

    def test_writes_the_report_with_the_old_file_name(self):
        self.assertTrue(self.call())
        self.assertTrue(os.path.exists(f"{self.tmp.name}/TIC 1_transits_validation_report.pdf"))

    def test_failure_is_swallowed_so_the_caller_can_fall_back(self):
        with mock.patch.object(lr, "write_landscape_report", side_effect=RuntimeError("boom")), \
                self.assertLogs(level=logging.ERROR):
            self.assertFalse(self.call())

    def test_system_exit_is_also_contained(self):
        with mock.patch.object(lr, "write_landscape_report", side_effect=SystemExit(1)), \
                self.assertLogs(level=logging.ERROR):
            self.assertFalse(self.call())


if __name__ == "__main__":
    unittest.main()
