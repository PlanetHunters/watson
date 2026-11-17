"""Characterization of the current (portrait) vetting report, before the landscape redesign."""
import glob
import os
import tempfile
import unittest

import pypdfium2 as pdfium

from watson.report import Report
from watson.tests.vetting_fixtures import make_vetting_dir, PERIOD, T0, DURATION_MIN, DEPTH_PPT


def pdf_info(path):
    doc = pdfium.PdfDocument(path)
    text = "\n".join(doc[i].get_textpage().get_text_range() for i in range(len(doc)))
    w, h = doc[0].get_size()
    return len(doc), text, (w, h)


class LegacyReportTest(unittest.TestCase):
    def build(self, **kw):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        t0s = make_vetting_dir(tmp.name, **{k: v for k, v in kw.items() if k != "summary"})
        summary = [0, 1] if kw.get("summary") else None
        Report(tmp.name, "report.pdf", "TIC 1", 63.37, -69.22, T0, PERIOD, DURATION_MIN, DEPTH_PPT, t0s, summary,
               12.8, 9.4, 8.8, 8.6, with_tpfs=False, is_summary=bool(summary)).create_report()
        return tmp.name, pdf_info(f"{tmp.name}/TIC 1_report.pdf")

    def test_full_report_is_portrait_with_known_sections(self):
        _, (pages, text, (w, h)) = self.build()
        self.assertLess(w, h)
        for section in ("Introduction", "Metrics summary", "WATSON-Net explainability", "TRICERATOPS results",
                        "Vetting information"):
            self.assertIn(section, text)
        self.assertGreaterEqual(pages, 6)

    def test_metrics_table_lists_every_metric_with_rounded_score(self):
        _, (_, text, _) = self.build()
        for metric in ("WATSON-NET", "Triceratops_FPP", "Bayesian_FPP", "snr", "bootstrap_fap"):
            self.assertIn(metric, text)

    def test_one_figure_per_single_transit_in_full_report(self):
        _, (_, text, _) = self.build(n_transits=4)
        self.assertEqual(sum(f"single transit no. {i}" in text for i in range(4)), 4)

    def test_summary_report_keeps_only_selected_transits(self):
        _, (_, text, _) = self.build(n_transits=4, summary=True)
        self.assertIn("single transit no. 0", text)
        self.assertNotIn("single transit no. 3", text)

    def test_missing_optional_inputs_do_not_break_the_report(self):
        _, (pages, text, _) = self.build(iatson=False, triceratops=False)
        self.assertNotIn("WATSON-Net explainability", text)
        self.assertNotIn("TRICERATOPS results", text)
        self.assertGreaterEqual(pages, 2)

    def test_missing_transit_image_is_skipped(self):
        tmp = tempfile.TemporaryDirectory(); self.addCleanup(tmp.cleanup)
        t0s = make_vetting_dir(tmp.name, n_transits=3)
        os.remove(glob.glob(f"{tmp.name}/single_transit_1_*.png")[0])
        Report(tmp.name, "r.pdf", "TIC 1", None, None, T0, PERIOD, DURATION_MIN, DEPTH_PPT, t0s, None,
               None, None, None, None, with_tpfs=False).create_report()
        _, text, _ = pdf_info(f"{tmp.name}/TIC 1_r.pdf")
        self.assertNotIn("single transit no. 1", text)
        self.assertIn("single transit no. 2", text)


if __name__ == "__main__":
    unittest.main()
