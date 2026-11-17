"""Synthetic vetting data directories shared by the report tests (no network, no TPFs needed)."""
import os

import numpy as np
import pandas as pd
from PIL import Image

PERIOD, T0, DURATION_MIN, DEPTH_PPT = 3.2899, 1327.51, 199.0, 6.082


def _png(path, size=(300, 200), color=(200, 200, 230)):
    Image.new("RGB", size, color).save(path)


def make_vetting_dir(directory, n_transits=4, iatson=True, triceratops=True, metrics=True, tpfs=False):
    """Create the files `Report.create_report` reads. Returns the list of transit t0s."""
    os.makedirs(directory, exist_ok=True)
    if metrics:
        pd.DataFrame({"metric": ["snr", "bootstrap_fap", "odd_even_depth_sigma", "secondary_snr"],
                      "score": [12.6, 0.0, 0.32, -0.47],
                      "passed": [True, True, np.nan, False]}).to_csv(f"{directory}/metrics.csv", index=False)
    if iatson:
        pd.DataFrame({"prediction_value_cal_mean": [0.994], "prediction_value_cal_std": [0.011]}) \
            .to_csv(f"{directory}/iatson_averages.csv", index=False)
        pd.DataFrame({"prediction_value_cal_mean": [0.9]}).to_csv(f"{directory}/iatson_predictions.csv", index=False)
        pd.DataFrame({"object_id": ["global_view", "odd_even"], "prediction_value_cal_mean": [0.12, -0.03]}) \
            .to_csv(f"{directory}/iatson_explain_branches.csv", index=False)
        pd.DataFrame({"object_id": ["period", "depth"], "prediction_value_cal_mean": [0.002, -0.001]}) \
            .to_csv(f"{directory}/iatson_explain_values.csv", index=False)
    if triceratops:
        tdir = f"{directory}/triceratops"
        os.makedirs(tdir, exist_ok=True)
        pd.DataFrame({"scenario": ["MEAN"], "FPP": [0.072], "FPP2": [0.01], "FPP3+": [0.005], "FPP_sys": [0.07],
                      "FPP2_sys": [0.01], "FPP3+_sys": [0.005], "NFPP": [2e-9]}).to_csv(f"{tdir}/validation.csv", index=False)
        pd.DataFrame({"combined_fpp": [5e-5], "combined_nfpp": [1e-9]}).to_csv(f"{tdir}/bayesian_fpp.csv", index=False)
        pd.DataFrame({"ID": [1], "scenario": ["TP"], "M_s": [0.5], "R_s": [0.5], "P_orb": [3.3], "inc": [88.0],
                      "b": [0.2], "ecc": [0.0], "w": [0.0], "R_p": [1.0], "M_EB": [0.0], "R_EB": [0.0], "prob": [0.9]}) \
            .to_csv(f"{tdir}/validation_scenarios.csv", index=False)
    _png(f"{directory}/odd_even_folded_curves.png", (400, 560))
    _png(f"{directory}/transit_depths.png", (560, 320))
    t0s = [round(T0 + i * PERIOD, 6) for i in range(n_transits)]
    for i, t in enumerate(t0s):
        _png(f"{directory}/single_transit_{i}_T0_{t}.png", (680, 960 if tpfs else 400))
    return t0s


def make_curve_and_transits(directory, n_transits=12, seed=1):
    """A noisy light curve with a box transit, plus the per-transit depth table. Returns (lc_file, transits)."""
    rng = np.random.default_rng(seed)
    time = np.arange(T0 - 2, T0 + (n_transits + 1) * PERIOD, 2 / 60 / 24)
    flux = 1 + rng.normal(0, 8e-4, time.size)
    phase = ((time - T0 + PERIOD / 2) % PERIOD) - PERIOD / 2
    flux[np.abs(phase) < DURATION_MIN / 60 / 24 / 2] -= DEPTH_PPT / 1000
    lc_file = f"{directory}/lc.csv"
    pd.DataFrame({"time": time, "flux": flux, "flux_err": 8e-4}).to_csv(lc_file, index=False)
    depths = DEPTH_PPT + rng.normal(0, 0.3, n_transits)
    depths[3] += 2.5  # a clear outlier
    transits = {"t0": [T0 + i * PERIOD for i in range(n_transits)], "depth": depths.tolist(),
                "depth_err": [0.3] * n_transits}
    return lc_file, transits


def make_landscape_dir(directory, n_transits=12, **kw):
    """A vetting directory plus the curve and transits the landscape report reads."""
    make_vetting_dir(directory, n_transits=n_transits, **kw)
    pd.DataFrame({"time": [T0, T0 + PERIOD], "time_folded": [.5, .5], "centroids_ra": [0.0, 0.1],
                  "centroids_dec": [0.0, 0.1]}).to_csv(f"{directory}/centroids.csv")
    lc_file, transits = make_curve_and_transits(directory, n_transits)
    pd.DataFrame({"time": [T0], "time_folded": [.5], "core_flux": [1.0], "halo_flux": [1.0], "og_flux": [0.0]}) \
        .to_csv(f"{directory}/og_dg.csv")
    pd.DataFrame({"name": ["diff_img", "mean"], "ra": [63.37, 63.37], "dec": [-69.22, -69.22], "ra_err": [.01, .01],
                  "dec_err": [.01, .01]}).to_csv(f"{directory}/source_offsets.csv")
    pd.DataFrame({"ID": [1, 2], "Tmag": [10.0, 12.0], "sep (arcsec)": [0.0, 30.0], "PA (E of N)": [0.0, 90.0]}) \
        .to_csv(f"{directory}/triceratops/stars.csv", index=False)
    return lc_file, transits


def pdf_text(path):
    """(pages, text, (width, height)) of a PDF. Skips the test when no PDF reader is installed."""
    import unittest
    try:
        import pypdfium2 as pdfium
    except ImportError:
        raise unittest.SkipTest("pypdfium2 is not installed")
    doc = pdfium.PdfDocument(path)
    text = "\n".join(doc[i].get_textpage().get_text_range() for i in range(len(doc)))
    return len(doc), text, doc[0].get_size()
