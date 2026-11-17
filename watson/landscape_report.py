"""Compact landscape vetting report.

Replaces the 80+ page portrait report by a handful of dense landscape pages:

1. verdict, the signal (folded, secondary, odd/even, depths) and the metrics,
2. WATSON-NET and TRICERATOPS,
3. is the signal on the target (source offset, centroids, optical ghost, per-sector pixels),
4. every single transit as a thumbnail,
5. detailed sheets (two per page) of the flagged transits only.

Everything is drawn with matplotlib ``Figure`` objects (no pyplot, safe to call from anywhere) and written with
``PdfPages``. Each page is built in isolation: a page that cannot be built becomes a short "unavailable" page
instead of breaking the report.
"""
import glob
import logging
import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle

INK, MUTED, GRID = "#1f2933", "#6b7785", "#e3e7ec"
OK, BAD, BLUE, ORANGE, RED, GOLD, INFO = "#2f8f5b", "#c0392b", "#2a6f97", "#e08a2b", "#d1495b", "#d9a21b", "#9aa5b1"
TRANSIT_BAND = "#f4c7c3"
PAGE_SIZE = (11.69, 8.27)
OUTLIER_SIGMA = 2.5
MAX_SHEETS = 12
THUMBNAILS_MAX = 240
TESS_PIXEL_ARCSEC = 21.0
PIXEL_MAPS_FRACTION = 0.30  # bottom part of folded_tpf_*.png: differential image and per-pixel BLS SNR maps
PLANET_SCENARIOS = ("TP", "PTP", "STP", "DTP", "BTP")
RC = {"font.size": 7, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
      "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 8, "axes.titleweight": "bold"}

METRIC_GROUPS = [
    ("Signal", [("snr", "SNR"), ("fast_snr", "SNR fast"), ("short_snr", "SNR short"), ("long_snr", "SNR long"),
                ("bootstrap_fap", "Bootstrap FAP")]),
    ("Plausibility", [("temp_stat", "Temperature"), ("albedo_stat", "Albedo")]),
    ("Secondary / odd-even", [("secondary_snr", "Secondary SNR"), ("odd_even_depth_diff", "Odd-even diff."),
                              ("odd_even_depth_sigma", "Odd-even σ"), ("odd_even_correlation", "Odd-even corr.")]),
    ("Centroids / offset", [("transit_offset_pos", "Offset"), ("transit_offset_err", "Offset err."),
                            ("centroids_ra_snr", "Centroid RA SNR"), ("centroids_dec_snr", "Centroid Dec SNR")]),
    ("Optical ghost", [("og_score", "Ghost score"), ("core_flux_snr", "Core SNR"), ("halo_flux_snr", "Halo SNR")]),
    ("Classifiers", [("WATSON-NET", "WATSON-NET"), ("WATSON-NET err", "WATSON-NET err."),
                     ("Triceratops_FPP", "TRIC. FPP"), ("Triceratops_NFPP", "TRIC. NFPP"),
                     ("Bayesian_FPP", "Combined FPP"), ("Bayesian_NFPP", "Combined NFPP")]),
]


# ---------------------------------------------------------------------------------------------- inputs
@dataclass
class ReportInputs:
    """Everything the report needs. Paths may be missing: the matching panels are skipped."""
    data_dir: str
    object_id: str
    ra: Optional[float]
    dec: Optional[float]
    t0: float
    period: float
    duration_min: float
    depth_ppt: float
    v: Optional[float] = None
    j: Optional[float] = None
    h: Optional[float] = None
    k: Optional[float] = None
    lc_file: Optional[str] = None
    lc_data_file: Optional[str] = None
    tpfs_dir: Optional[str] = None
    apertures_file: Optional[str] = None
    star_file: Optional[str] = None
    transits: Optional[dict] = None
    transits_mask: Optional[list] = field(default=None)

    @property
    def duration_h(self):
        return self.duration_min / 60


def _read_csv(path, **kw):
    return pd.read_csv(path, **kw) if path and os.path.exists(path) else None


# ---------------------------------------------------------------------------------------------- numerics
def binned(x, y, lo, hi, n, robust=False, minimum=5):
    """Mean (or median) and standard error of ``y`` in ``n`` equal bins of ``x`` between ``lo`` and ``hi``."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    edges = np.linspace(lo, hi, n + 1)
    idx = np.digitize(x, edges) - 1
    cx, cy, ce = [], [], []
    for i in range(n):
        s = (idx == i) & np.isfinite(y)
        if s.sum() >= minimum:
            cx.append(x[s].mean())
            cy.append(np.median(y[s]) if robust else y[s].mean())
            ce.append(y[s].std() / np.sqrt(s.sum()))
    return np.array(cx), np.array(cy), np.array(ce)


def box_model(x_hours, depth, duration_h, ingress_h=0.15):
    """Trapezoid transit model on a time axis in hours from the mid-transit."""
    x = np.abs(np.asarray(x_hours, dtype=float))
    ingress_h = min(ingress_h, duration_h / 2)
    y = np.ones_like(x)
    y[x < duration_h / 2 - ingress_h] = 1 - depth
    edge = (x >= duration_h / 2 - ingress_h) & (x < duration_h / 2)
    y[edge] = 1 - depth * (duration_h / 2 - x[edge]) / ingress_h
    return y


def phase_hours(time, t0, period, shift=0.0):
    """Hours from the closest mid-transit (``shift`` = fraction of the period, 0.5 for the secondary)."""
    return (((np.asarray(time) - t0 - shift * period + period / 2) % period) - period / 2) * 24


def transit_depths_summary(transits, depth_ppt):
    """Per-transit depth statistics: z-scores against the weighted mean and the flagged transit indexes.

    Returns ``(frame, z, flagged)`` where ``flagged`` maps the index of a transit in ``frame`` to the reason.
    Flagged transits are the deepest, the shallowest, the closest to the candidate depth and every
    ``|depth - mean| > 2.5 sigma`` outlier (the rule of the old summary report plus the outliers).
    """
    frame = pd.DataFrame({"t0": np.asarray(transits["t0"], float), "depth": np.asarray(transits["depth"], float),
                          "depth_err": np.asarray(transits["depth_err"], float)})
    frame = frame[np.isfinite(frame.depth) & np.isfinite(frame.depth_err) & (frame.depth_err > 0)].reset_index(drop=True)
    if len(frame) == 0:
        return frame, np.array([]), {}
    mean = np.average(frame.depth, weights=1 / frame.depth_err ** 2)
    z = ((frame.depth - mean) / frame.depth_err).to_numpy()
    reasons = {}
    if len(frame) > 2:
        reasons.setdefault(int(frame.depth.idxmax()), []).append("deepest")
        reasons.setdefault(int(frame.depth.idxmin()), []).append("shallowest")
        reasons.setdefault(int((frame.depth - depth_ppt).abs().idxmin()), []).append("closest to the candidate depth")
    for i in np.where(np.abs(z) > OUTLIER_SIGMA)[0]:
        reasons.setdefault(int(i), []).append(f"depth outlier ({z[i]:+.1f} σ)")
    flagged = {i: ", ".join(r) for i, r in sorted(reasons.items())}
    return frame, z, flagged


# ---------------------------------------------------------------------------------------------- metrics
def _row(metric, score, passed):
    return {"metric": metric, "score": score, "passed": passed}


def classifier_metrics(data_dir):
    """Metrics of WATSON-NET, TRICERATOPS and the combined FPP with the same verdict rules as the old report."""
    rows = []
    averages = _read_csv(f"{data_dir}/iatson_averages.csv")
    if averages is not None and len(averages):
        mean, std = float(averages.loc[0, "prediction_value_cal_mean"]), float(averages.loc[0, "prediction_value_cal_std"])
        rows.append(_row("WATSON-NET", mean, True if mean + std >= 0.646 else (False if mean + std <= 0.056 else np.nan)))
        rows.append(_row("WATSON-NET err", std, True if std < 0.015 else (np.nan if std < 0.1 else False)))
    validation = _read_csv(f"{data_dir}/triceratops/validation.csv")
    if validation is not None and (validation.scenario == "MEAN").any():
        mean_row = validation[validation.scenario == "MEAN"].iloc[0]
        fpp, nfpp = float(mean_row.FPP), float(mean_row.NFPP)
        rows.append(_row("Triceratops_FPP", fpp, True if fpp <= 0.015 else (np.nan if fpp <= 0.5 else False)))
        rows.append(_row("Triceratops_NFPP", nfpp, True if nfpp <= 0.001 else (np.nan if nfpp <= 0.1 else False)))
    bayes = _read_csv(f"{data_dir}/triceratops/bayesian_fpp.csv")
    if bayes is not None and len(bayes):
        fpp, nfpp = float(bayes.loc[0, "combined_fpp"]), float(bayes.loc[0, "combined_nfpp"])
        rows.append(_row("Bayesian_FPP", fpp, True if fpp <= 0.01 else (np.nan if fpp <= 0.5 else False)))
        rows.append(_row("Bayesian_NFPP", nfpp, True if nfpp <= 0.001 else (np.nan if nfpp <= 0.1 else False)))
    return pd.DataFrame(rows, columns=["metric", "score", "passed"])


def load_metrics(data_dir):
    """All the metrics (``metrics.csv`` plus the classifiers) with ``passed`` as 1 / 0 / NaN (no verdict)."""
    frames = [f for f in (_read_csv(f"{data_dir}/metrics.csv"), classifier_metrics(data_dir))
              if f is not None and len(f)]
    if not frames:
        return pd.DataFrame(columns=["metric", "score", "passed"]).set_index("metric")
    metrics = pd.concat([f[["metric", "score", "passed"]] for f in frames], ignore_index=True)
    metrics["passed"] = pd.to_numeric(metrics["passed"].replace({"True": 1, "False": 0, True: 1, False: 0}),
                                      errors="coerce")
    metrics["score"] = pd.to_numeric(metrics["score"], errors="coerce")
    return metrics.drop_duplicates("metric", keep="first").set_index("metric")


def status_of(metrics, name):
    """'pass', 'fail' or 'none' (no verdict, or the metric is missing)."""
    if name not in metrics.index or pd.isna(metrics.loc[name, "passed"]):
        return "none"
    return "pass" if metrics.loc[name, "passed"] == 1 else "fail"


def verdict(metrics):
    """Headline, colour and detail of the report: it states what failed and never claims to confirm a planet."""
    failed = [m for m in metrics.index if status_of(metrics, m) == "fail"]
    passed = [m for m in metrics.index if status_of(metrics, m) == "pass"]
    undecided = [m for m in metrics.index if status_of(metrics, m) == "none"]
    detail = f"{len(passed)} of {len(passed) + len(failed)} tests with a verdict passed · {len(undecided)} without a verdict"
    if failed:
        return f"{len(failed)} test{'s' if len(failed) > 1 else ''} failed", BAD, detail + " · failed: " + ", ".join(failed)
    if passed:
        return "No test failed", OK, detail
    return "No verdict available", INFO, detail


def _fmt(value):
    if value is None or not np.isfinite(value):
        return "-"
    if value != 0 and (abs(value) < 1e-3 or abs(value) >= 1e4):
        return f"{value:.1e}"
    return f"{value:.3g}"


# ---------------------------------------------------------------------------------------------- drawing helpers
def _axes(fig, x, y, w, h, **kw):
    return fig.add_axes([x, y, w, h], **kw)


def _title(ax, text, note=None):
    ax.set_title(text, loc="left", pad=3)
    if note:
        ax.text(1, 1.025, note, transform=ax.transAxes, ha="right", fontsize=6.5, color=MUTED)


def _head(fig, inputs, title, subtitle, page):
    fig.text(.03, .962, title, fontsize=12.5, fontweight="bold", color=INK)
    fig.text(.03, .940, subtitle, fontsize=7.5, color=MUTED)
    fig.text(.97, .962, f"{inputs.object_id} · page {page}", fontsize=7, color=MUTED, ha="right")
    fig.add_artist(Line2D([.03, .97], [.93, .93], color=GRID, lw=.8))


def _unavailable(ax, text="not available"):
    ax.axis("off")
    ax.text(.5, .5, text, ha="center", va="center", color=MUTED, fontsize=7)


def _new_figure():
    return Figure(figsize=PAGE_SIZE)


def _band(ax, half_width):
    ax.axvspan(-half_width, half_width, color=TRANSIT_BAND, alpha=.5, lw=0)


def _star_line(inputs):
    parts = [f"P = {inputs.period:.5f} d", f"T0 = {inputs.t0:.4f}", f"duration {inputs.duration_h:.2f} h",
             f"depth {inputs.depth_ppt:.3f} ppt"]
    if inputs.ra is not None and inputs.dec is not None:
        from astropy import units as u
        from astropy.coordinates import Angle
        parts.append("RA " + Angle(inputs.ra, u.deg).to_string(unit=u.hourangle, sep=":", precision=1) +
                     "  Dec " + Angle(inputs.dec, u.deg).to_string(unit=u.deg, sep=":", precision=1))
    mags = [f"{name} {value:.1f}" for name, value in zip("VJHK", (inputs.v, inputs.j, inputs.h, inputs.k))
            if value is not None]
    return " · ".join(parts + ([" ".join(mags)] if mags else []))


def _lc(inputs):
    frame = _read_csv(inputs.lc_file)
    if frame is None:
        return None
    frame = frame.sort_values("time")
    return frame.time.to_numpy(float), frame.flux.to_numpy(float)


# ---------------------------------------------------------------------------------------------- page 1
def page_verdict(inputs, metrics):
    """Verdict, folded transit, secondary, odd/even, single-transit depths and every metric."""
    fig = _new_figure()
    _head(fig, inputs, "Verdict and signal", _star_line(inputs), 1)
    headline, colour, detail = verdict(metrics)
    ax = _axes(fig, .03, .845, .94, .075)
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=.01",
                                fc="#e6f3ec" if colour == OK else ("#fbe7e4" if colour == BAD else "#eef1f4"), ec="none"))
    ax.text(.012, .66, headline, fontsize=15, fontweight="bold", color=colour, va="center")
    ax.text(.012, .22, detail if len(detail) < 150 else detail[:147] + "...", fontsize=7.5, va="center")
    tiles = []
    if "WATSON-NET" in metrics.index:
        tiles.append(("WATSON-NET", f"{metrics.loc['WATSON-NET', 'score']:.3f}", "planet probability", "WATSON-NET"))
    if "Triceratops_FPP" in metrics.index:
        tiles.append(("TRICERATOPS FPP", f"{metrics.loc['Triceratops_FPP', 'score'] * 100:.1f} %", "", "Triceratops_FPP"))
    if "Bayesian_FPP" in metrics.index:
        tiles.append(("Combined FPP", _fmt(metrics.loc["Bayesian_FPP", "score"]), "NN + TRICERATOPS", "Bayesian_FPP"))
    if "transit_offset_pos" in metrics.index:
        tiles.append(("Source offset", "OK" if status_of(metrics, "transit_offset_pos") == "pass" else "check",
                      "from target", "transit_offset_pos"))
    tile_colour = {"pass": OK, "fail": BAD, "none": INFO}
    for i, (label, big, small, key) in enumerate(tiles):
        x = .43 + i * .14
        ax.text(x, .76, label, fontsize=6.5, color=MUTED, va="center")
        ax.text(x, .46, big, fontsize=14, fontweight="bold", color=tile_colour[status_of(metrics, key)], va="center")
        ax.text(x, .15, small, fontsize=6.5, color=MUTED, va="center")
    curve = _lc(inputs)
    half = max(3.4, 2.5 * inputs.duration_h)
    depth = inputs.depth_ppt / 1000
    if curve is not None:
        t, f = curve
        ph = phase_hours(t, inputs.t0, inputs.period)
        sec = phase_hours(t, inputs.t0, inputs.period, 0.5)
        epoch = np.round((t - inputs.t0) / inputs.period).astype(int)
        sel = np.abs(ph) < half
        ylim = (1 - 2.2 * depth, 1 + 1.4 * depth)
        a, ar = _axes(fig, .06, .595, .20, .20), None
        ar = _axes(fig, .06, .505, .20, .07, sharex=a)
        a.scatter(ph[sel][::5], f[sel][::5], s=.8, color=MUTED, alpha=.1, lw=0, rasterized=True)
        bx, by, be = binned(ph[sel], f[sel], -half, half, 36)
        a.errorbar(bx, by, be, fmt="o", ms=2.2, color=INK, lw=.6)
        xm = np.linspace(-half, half, 300)
        a.plot(xm, box_model(xm, depth, inputs.duration_h), color=BAD, lw=1.3)
        a.set_ylim(*ylim)
        a.grid(color=GRID, lw=.4)
        a.tick_params(labelbottom=False)
        a.set_ylabel("flux")
        _title(a, "Folded transit")
        ar.errorbar(bx, (by - box_model(bx, depth, inputs.duration_h)) * 1e3, be * 1e3, fmt="o", ms=2, color=INK, lw=.5)
        ar.axhline(0, color=BAD, lw=.7)
        ar.set_ylabel("res. ppt", fontsize=6)
        ar.set_xlabel("hours from mid-transit")
        ar.grid(color=GRID, lw=.4)
        a = _axes(fig, .325, .505, .215, .29)
        sel2 = np.abs(sec) < 2 * half
        a.scatter(sec[sel2][::5], f[sel2][::5], s=.8, color=MUTED, alpha=.1, lw=0, rasterized=True)
        b2 = binned(sec[sel2], f[sel2], -2 * half, 2 * half, 32)
        a.errorbar(b2[0], b2[1], b2[2], fmt="o", ms=2.2, color=INK, lw=.6)
        _band(a, inputs.duration_h / 2)
        a.axhline(1, color=BAD, lw=.7, ls="--")
        a.set_ylim(*ylim)
        a.set_xlabel("hours from phase 0.5")
        a.grid(color=GRID, lw=.4)
        sec_snr = metrics.loc["secondary_snr", "score"] if "secondary_snr" in metrics.index else None
        _title(a, "Secondary event", None if sec_snr is None else f"SNR {sec_snr:+.2f}")
        a = _axes(fig, .605, .505, .175, .29)
        for parity, color, label in ((1, RED, "odd"), (0, BLUE, "even")):
            s = sel & (epoch % 2 == parity)
            b3 = binned(ph[s], f[s], -half, half, 22)
            a.errorbar(b3[0], b3[1], b3[2], fmt="o-", ms=2.2, lw=.6, color=color, label=label)
        a.plot(xm, box_model(xm, depth, inputs.duration_h), color=INK, lw=1, ls=":")
        a.set_ylim(*ylim)
        a.set_xlabel("hours from mid-transit")
        a.grid(color=GRID, lw=.4)
        a.legend(frameon=False, fontsize=6, loc="upper right", ncol=2)
        note = None
        if "odd_even_depth_sigma" in metrics.index:
            note = f"{metrics.loc['odd_even_depth_sigma', 'score']:.2f} σ"
        _title(a, "Odd / even folded", note)
    else:
        for rect in ((.06, .505, .20, .29), (.325, .505, .215, .29), (.605, .505, .175, .29)):
            _unavailable(_axes(fig, *rect), "light curve not available")
    a = _axes(fig, .835, .505, .135, .29)
    if inputs.transits is not None:
        frame, z, _ = transit_depths_summary(inputs.transits, inputs.depth_ppt)
        epoch_of = np.round((frame.t0.to_numpy() - inputs.t0) / inputs.period).astype(int)
        w, d = (1 / frame.depth_err ** 2).to_numpy(), frame.depth.to_numpy()

        def weighted(mask):
            return (float(np.sum(w[mask] * d[mask]) / np.sum(w[mask])), float(1 / np.sqrt(np.sum(w[mask])))) \
                if mask.any() else (np.nan, np.nan)

        groups = [weighted(np.ones(len(d), bool)), weighted(epoch_of % 2 == 1), weighted(epoch_of % 2 == 0)]
        a.errorbar(range(3), [g[0] for g in groups], [g[1] for g in groups], fmt="none", ecolor=INK, capsize=3)
        for i, c in enumerate((INK, RED, BLUE)):
            a.plot(i, groups[i][0], "o", ms=5, color=c)
        a.set_xticks(range(3))
        a.set_xticklabels(["all", "odd", "even"])
        a.set_xlim(-.5, 2.5)
        a.set_ylabel("depth (ppt)")
        a.grid(axis="y", color=GRID, lw=.4)
        note = None
        if np.isfinite(groups[1][0]) and np.isfinite(groups[2][0]):
            note = f"Δ {groups[1][0] - groups[2][0]:.2f} ± {np.hypot(groups[1][1], groups[2][1]):.2f}"
        _title(a, "Odd vs even depth", note)
    else:
        _unavailable(a, "single transits not available")
    _metrics_columns(fig, metrics)
    _page1_bottom(fig, inputs)
    fig.text(.03, .012, "● passed   ● failed   ● no verdict.  Verdicts come from the thresholds of each test in WATSON "
                        "(SNR > 3, |centroid SNR| < 3, odd-even < 1σ passes and > 3σ fails, FPP ≤ 0.015, NFPP ≤ 0.001...).",
             fontsize=6.5, color=MUTED)
    return fig


def _metrics_columns(fig, metrics):
    dot = {"pass": OK, "fail": BAD, "none": INFO}
    columns = []
    for name, items in METRIC_GROUPS:
        present = [(k, label) for k, label in items if k in metrics.index]
        if present:
            columns.append((name, present))
    if not columns:
        return
    width = .94 / len(columns)
    for gi, (name, items) in enumerate(columns):
        x, y = .035 + gi * width, .425
        fig.text(x, y, name.upper(), fontsize=6.5, color=MUTED, fontweight="bold")
        fig.add_artist(Line2D([x, x + width - .02], [y - .008, y - .008], color=GRID, lw=.6))
        y -= .03
        for key, label in items:
            fig.text(x, y, label, fontsize=7.2)
            fig.text(x + width - .045, y, _fmt(metrics.loc[key, "score"]), fontsize=7.2, fontweight="bold", ha="right")
            fig.patches.append(Circle((x + width - .03, y + .005), .0042, color=dot[status_of(metrics, key)],
                                      transform=fig.transFigure))
            y -= .024


def _page1_bottom(fig, inputs):
    """Candidate and star table, plus the FPP/NFPP map when TRICERATOPS ran."""
    star = _read_csv(inputs.star_file)
    fig.text(.605, .245, "CANDIDATE AND STAR", fontsize=6.5, color=MUTED, fontweight="bold")
    rows = [("Period", f"{inputs.period:.6f} d"), ("T0", f"{inputs.t0:.4f}"), ("Duration", f"{inputs.duration_h:.2f} h"),
            ("Depth", f"{inputs.depth_ppt:.3f} ppt")]
    if star is not None and len(star):
        s = star.iloc[0]
        parts = []
        for col, label, fmt in (("Teff_star", "Teff", "{:.0f} K"), ("R_star", "R", "{:.3f} R☉"), ("M_star", "M", "{:.3f} M☉")):
            if col in s and np.isfinite(s[col]):
                parts.append(f"{label} {fmt.format(s[col])}")
        if parts:
            rows.append(("Star", " · ".join(parts)))
    for i, (k, v) in enumerate(rows):
        fig.text(.605, .215 - i * .028, k, fontsize=7.5, color=MUTED)
        fig.text(.69, .215 - i * .028, v, fontsize=7.5, fontweight="bold")
    validation = _read_csv(f"{inputs.data_dir}/triceratops/validation.csv")
    if validation is not None:
        _fpp_map(_axes(fig, .055, .075, .205, .175), validation)
    stars = _read_csv(f"{inputs.data_dir}/triceratops/stars.csv")
    if stars is not None and len(stars) > 1:
        _field_map(_axes(fig, .335, .075, .205, .175), stars, small=True)


def _fpp_map(a, validation):
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(1e-4, 1)
    a.set_ylim(1e-12, 1)
    a.add_patch(Rectangle((1e-4, 1e-12), .5, 1e-3 - 1e-12, color="#d6efdf", lw=0))
    a.text(2e-4, 3e-11, "likely planet", fontsize=6.5, color=OK)
    scen = validation[validation.scenario != "MEAN"]
    a.scatter(scen.FPP.clip(lower=1e-4), scen.NFPP.clip(lower=1e-12), s=14, color=BLUE, label="scenarios")
    mean = validation[validation.scenario == "MEAN"]
    if len(mean):
        a.scatter(mean.FPP.clip(lower=1e-4), mean.NFPP.clip(lower=1e-12), s=36, color=INK, marker="D", label="mean")
    a.legend(frameon=False, fontsize=6, loc="upper right")
    a.set_xlabel("FPP")
    a.set_ylabel("NFPP")
    a.grid(color=GRID, lw=.4)
    _title(a, "Validation map")


def _field_map(a, stars, small=False):
    a.set_aspect("equal")
    sep = stars["sep (arcsec)"].to_numpy(float)
    pa = np.radians(stars["PA (E of N)"].to_numpy(float))
    tmag = stars["Tmag"].to_numpy(float)
    size = np.clip((400 if small else 500) * 10 ** (-.4 * (tmag[1:] - tmag[0])), 10, 220)
    sc = a.scatter((-sep * np.sin(pa))[1:], (sep * np.cos(pa))[1:], s=size, c=tmag[1:], cmap="viridis_r", alpha=.85, lw=0)
    a.scatter([0], [0], s=80, marker="*", color=GOLD, zorder=3, ec=INK)
    reach = float(max(10, np.ceil(sep.max() / 10) * 10))
    a.set_xlim(reach, -reach)
    a.set_ylim(-reach, reach)
    a.set_xlabel("east (″)")
    a.set_ylabel("north (″)")
    a.grid(color=GRID, lw=.4)
    _title(a, f"Nearby stars ({len(stars) - 1})", None if small else "colour: Tmag")
    return sc


# ---------------------------------------------------------------------------------------------- page 2
def page_classifiers(inputs, metrics):
    """WATSON-NET branches and inputs, TRICERATOPS map, scenarios, field and table."""
    fig = _new_figure()
    _head(fig, inputs, "Classifiers · WATSON-NET and TRICERATOPS",
          "Neural network probability and how it is built · statistical validation against false-positive scenarios", 2)
    averages = _read_csv(f"{inputs.data_dir}/iatson_averages.csv")
    branches = _read_csv(f"{inputs.data_dir}/iatson_explain_branches.csv")
    values = _read_csv(f"{inputs.data_dir}/iatson_explain_values.csv")
    prefix = re.escape(str(inputs.object_id)) + r"\s*"
    if averages is not None and len(averages):
        mean, std = float(averages.loc[0, "prediction_value_cal_mean"]), float(averages.loc[0, "prediction_value_cal_std"])
        colour = {"pass": OK, "fail": BAD, "none": GOLD}[status_of(metrics, "WATSON-NET")]
        a = _axes(fig, .03, .74, .20, .16)
        a.axis("off")
        a.text(0, .95, "WATSON-NET planet probability", fontsize=7, color=MUTED, va="top")
        a.text(0, .45, f"{mean:.3f}", fontsize=28, fontweight="bold", color=colour, va="center")
        a.text(.55, .38, f"± {std:.3f}\ncalibrated", fontsize=7, color=MUTED, va="center")
        a.add_patch(Rectangle((0, .05), 1, .1, color=GRID, transform=a.transAxes))
        a.add_patch(Rectangle((0, .05), min(max(mean, 0), 1), .1, color=colour, transform=a.transAxes))
    if branches is not None and len(branches):
        b = branches.copy()
        b["name"] = b.object_id.astype(str).str.replace(prefix, "", regex=True).str.replace(" branch", "", regex=False)
        b = b.sort_values("prediction_value_cal_mean")
        a = _axes(fig, .115, .08, .13, .60)
        a.barh(range(len(b)), b.prediction_value_cal_mean, xerr=b.prediction_value_cal_std,
               color=[OK if v >= 0 else BAD for v in b.prediction_value_cal_mean], error_kw=dict(lw=.4, ecolor=MUTED), height=.7)
        a.set_yticks(range(len(b)))
        a.set_yticklabels(b["name"], fontsize=6)
        a.axvline(0, color=INK, lw=.6)
        a.grid(axis="x", color=GRID, lw=.4)
        a.set_xlabel("impact (positive: supports a planet)", fontsize=6.5)
        _title(a, "Network branches")
    if values is not None and len(values):
        v = values.copy()
        v["feature"] = v.object_id.astype(str).str.replace(prefix, "", regex=True).str.replace(r" value=.*", "", regex=True)
        v["setting"] = v.object_id.astype(str).str.extract(r"value=([\d.e+-]+)")[0].astype(float).round(1)
        piv = v.pivot_table(index="feature", columns="setting", values="prediction_value_cal_mean")
        if piv.shape[1] and piv.notna().any().any():
            piv = piv.loc[piv.abs().max(axis=1).sort_values().index]
            a = _axes(fig, .335, .08, .135, .60)
            y = np.arange(len(piv))
            n = piv.shape[1]
            for k, (name, color) in enumerate(zip(piv.columns, (BLUE, ORANGE, BAD, GOLD))):
                a.barh(y + (k - (n - 1) / 2) * (.78 / n), piv[name], height=.78 / n, color=color, label=f"set to {name:g}")
            a.axvline(0, color=INK, lw=.6)
            a.set_yticks(y)
            a.set_yticklabels(piv.index, fontsize=6)
            a.legend(frameon=False, fontsize=6, loc="lower left")
            a.grid(axis="x", color=GRID, lw=.4)
            a.set_xlabel("impact of changing the input", fontsize=6.5)
            _title(a, "Inputs")
    if averages is None:
        _unavailable(_axes(fig, .03, .08, .44, .82), "WATSON-NET was not run")
    validation = _read_csv(f"{inputs.data_dir}/triceratops/validation.csv")
    bayes = _read_csv(f"{inputs.data_dir}/triceratops/bayesian_fpp.csv")
    if validation is None:
        _unavailable(_axes(fig, .53, .08, .44, .82), "TRICERATOPS was not run")
        return fig
    mean_row = validation[validation.scenario == "MEAN"]
    headline = "TRICERATOPS"
    fig.text(.53, .915, headline, fontsize=8, fontweight="bold")
    if len(mean_row):
        text = f"FPP {mean_row.iloc[0].FPP * 100:.1f} %  ·  NFPP {mean_row.iloc[0].NFPP:.0e}"
        if bayes is not None and len(bayes):
            text += f"  ·  combined FPP {bayes.loc[0, 'combined_fpp']:.0e}"
        fig.text(.53, .895, text, fontsize=8, fontweight="bold",
                 color={"pass": OK, "fail": BAD, "none": GOLD}[status_of(metrics, "Triceratops_FPP")])
    _fpp_map(_axes(fig, .53, .52, .20, .34), validation)
    scen = _read_csv(f"{inputs.data_dir}/triceratops/validation_scenarios.csv")
    a = _axes(fig, .79, .52, .18, .34)
    if scen is not None and len(scen):
        top = scen.sort_values("prob", ascending=False).head(8).iloc[::-1]
        a.barh(top.scenario.astype(str), top.prob, color=[OK if s in PLANET_SCENARIOS else BAD for s in top.scenario], height=.7)
        a.set_xlabel("probability")
        a.grid(axis="x", color=GRID, lw=.4)
        _title(a, "Scenario probabilities", "green: planet")
    else:
        _unavailable(a)
    stars = _read_csv(f"{inputs.data_dir}/triceratops/stars.csv")
    a = _axes(fig, .53, .08, .17, .36)
    if stars is not None and len(stars) > 1:
        _field_map(a, stars)
    else:
        _unavailable(a, "no nearby stars")
    a = _axes(fig, .73, .08, .24, .36)
    a.axis("off")
    for x, text, ha in ((0, "SCENARIO", "left"), (.35, "FPP", "right"), (.62, "FPP2", "right"), (.92, "NFPP", "right")):
        a.text(x, 1, text, fontsize=6.5, color=MUTED, fontweight="bold", va="top", ha=ha)
    for i, (_, r) in enumerate(validation.head(9).iterrows()):
        y, bold = .88 - i * .1, ("bold" if r.scenario == "MEAN" else "normal")
        a.text(0, y, str(r.scenario), fontsize=7.5, va="top", fontweight=bold)
        a.text(.35, y, f"{r.FPP:.3f}", fontsize=7.5, va="top", ha="right", fontweight=bold)
        a.text(.62, y, f"{r.FPP2:.4f}", fontsize=7.5, va="top", ha="right", fontweight=bold)
        a.text(.92, y, f"{r.NFPP:.0e}", fontsize=7.5, va="top", ha="right", fontweight=bold)
    fig.text(.03, .018, "WATSON-NET: impact of each branch and input on the prediction. TRICERATOPS best-fit plots of "
                        "each scenario stay in the triceratops/ folder.", fontsize=6.5, color=MUTED)
    return fig


# ---------------------------------------------------------------------------------------------- page 3
def _trim_white(image):
    from PIL import ImageChops, Image
    image = image.convert("RGB")
    bbox = ImageChops.difference(image, Image.new("RGB", image.size, "white")).getbbox()
    return image.crop(bbox) if bbox else image


def _folded_phase_panel(ax, frame, column, label, color, ylabel, note, inputs, show_xlabel):
    if frame is None or column not in frame:
        _unavailable(ax)
        return
    ph = ((frame.time.to_numpy(float) - inputs.t0) / inputs.period + .5) % 1
    y = frame[column].to_numpy(float)
    sel = (ph > .4) & (ph < .6) & np.isfinite(y)
    if sel.sum() < 10:
        _unavailable(ax)
        return
    x, y = ph[sel], y[sel]
    ax.scatter(x[::4], y[::4], s=.8, color=MUTED, alpha=.1, lw=0, rasterized=True)
    bx, by, _ = binned(x, y, .4, .6, 40, robust=True, minimum=15)
    ax.plot(bx, by, "-o", ms=2, lw=1, color=color)
    half = inputs.duration_h / 24 / inputs.period / 2
    ax.axvspan(.5 - half, .5 + half, color=TRANSIT_BAND, alpha=.5, lw=0)
    ax.set_ylim(*np.percentile(y, [2, 98]))
    ax.set_xlabel("phase" if show_xlabel else "")
    ax.set_ylabel(ylabel, fontsize=6.5)
    ax.grid(color=GRID, lw=.4)
    _title(ax, label, note)


def _source_offset_panel(ax, inputs, tpf, aperture_pixels):
    from astropy.coordinates import SkyCoord
    img = np.nanmedian(tpf.flux.value, axis=0)
    c0, r0 = tpf.column, tpf.row
    ny, nx = img.shape
    ax.imshow(img, origin="lower", extent=[c0 - .5, c0 + nx - .5, r0 - .5, r0 + ny - .5], cmap="viridis")
    for c, r in aperture_pixels:
        ax.add_patch(Rectangle((c - .5, r - .5), 1, 1, fill=False, hatch="///", ec="#ffd166", lw=.7))

    def px(ra, dec):
        x, y = tpf.wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
        return c0 + float(x), r0 + float(y)

    if inputs.ra is not None:
        tx, ty = px(inputs.ra, inputs.dec)
        ax.plot(tx, ty, "*", ms=12, color="#ffd166", mec=INK, label="target", zorder=5)
    offsets = _read_csv(f"{inputs.data_dir}/source_offsets.csv")
    note = None
    if offsets is not None:
        offsets = offsets.set_index("name")
        for name, (marker, col, label) in {"diff_img": ("D", "#00bcd4", "differential image"),
                                           "px_bls": ("s", "#7bd88f", "per-pixel BLS"),
                                           "mean": ("o", RED, "mean offset")}.items():
            if name not in offsets.index:
                continue
            r = offsets.loc[name]
            x, y = px(r.ra, r.dec)
            ax.plot(x, y, marker, ms=6, color=col, mec="white", label=label, zorder=6)
            cosd = np.cos(np.radians(r.dec))
            ax.add_patch(Ellipse((x, y), 2 * r.ra_err * cosd * 3600 / TESS_PIXEL_ARCSEC, 2 * r.dec_err * 3600 / TESS_PIXEL_ARCSEC,
                                 fill=False, ec=col, lw=1, zorder=6))
        if "mean" in offsets.index and inputs.ra is not None:
            r = offsets.loc["mean"]
            dist = np.hypot((r.ra - inputs.ra) * np.cos(np.radians(inputs.dec)), r.dec - inputs.dec) * 3600
            note = f"{dist:.1f}″ ± {r.ra_err * 3600:.1f}″"
    ax.set_xlim(c0 + 1.5, c0 + nx - 2.5)
    ax.set_ylim(r0 + 1.5, r0 + ny - 2.5)
    ax.legend(fontsize=6, loc="upper left", framealpha=.85)
    ax.set_xlabel("pixel column")
    ax.set_ylabel("pixel row")
    _title(ax, f"Source offset (sector {_sector_of(tpf)})", note)


def _sector_of(tpf):
    for attribute in ("sector", "campaign", "quarter"):
        value = getattr(tpf, attribute, None)
        if value is not None:
            return value
    return "?"


def _aperture_pixels(apertures, tpf):
    from watson.watson import Watson
    if not isinstance(apertures, dict):
        return []
    author, cadence = Watson.get_author_cadence_from_tpf_name(tpf)
    return [tuple(p) for p in Watson.get_aperture_for_sector(apertures, _sector_of(tpf), author, cadence)]


def _load_apertures(inputs):
    if inputs.apertures_file and os.path.exists(inputs.apertures_file):
        import yaml
        with open(inputs.apertures_file) as handle:
            return yaml.load(handle, yaml.SafeLoader)["sectors"]
    return None


def _load_tpfs(inputs):
    from watson.watson import Watson
    if not (inputs.lc_file and inputs.tpfs_dir and os.path.exists(inputs.tpfs_dir)):
        return None, None, [], None
    lc, lc_data, lc_data_norm, tpfs = Watson.initialize_lc_and_tpfs(
        inputs.object_id, inputs.lc_file, inputs.lc_data_file or "", inputs.tpfs_dir, transits_mask=inputs.transits_mask)
    return lc, lc_data, tpfs, lc_data_norm


def page_on_target(inputs, metrics):
    """Source offset, centroid shifts, optical ghost and the per-sector pixel analysis."""
    fig = _new_figure()
    _head(fig, inputs, "Is the signal on the target?",
          "Source offset, centroid shifts and optical ghost around the transit, and the pixel analysis of every sector", 3)
    _, _, tpfs, _ = _load_tpfs(inputs)
    apertures = _load_apertures(inputs)
    a = _axes(fig, .045, .46, .27, .43)
    if tpfs:
        _source_offset_panel(a, inputs, tpfs[0], _aperture_pixels(apertures, tpfs[0]))
    else:
        _unavailable(a, "no target pixel files")
    cen = _read_csv(f"{inputs.data_dir}/centroids.csv")
    og = _read_csv(f"{inputs.data_dir}/og_dg.csv")

    def note(name, fmt="{:+.2f}"):
        return None if name not in metrics.index else "SNR " + fmt.format(metrics.loc[name, "score"])

    _folded_phase_panel(_axes(fig, .37, .66, .27, .23), cen, "centroids_ra", "Centroid shift RA", BLUE, "drift",
                        note("centroids_ra_snr"), inputs, False)
    _folded_phase_panel(_axes(fig, .70, .66, .27, .23), cen, "centroids_dec", "Centroid shift Dec", BLUE, "drift",
                        note("centroids_dec_snr"), inputs, False)
    _folded_phase_panel(_axes(fig, .37, .46, .27, .15), og, "core_flux", "Ghost: core flux", ORANGE, "flux",
                        note("core_flux_snr", "{:.1f}"), inputs, True)
    _folded_phase_panel(_axes(fig, .70, .46, .27, .15), og, "halo_flux", "Ghost: halo flux", ORANGE, "flux",
                        note("halo_flux_snr", "{:.2f}"), inputs, True)
    images = sorted(glob.glob(f"{inputs.data_dir}/folded_tpf_*.png"))
    if images:
        from PIL import Image
        fig.text(.03, .385, "PIXEL ANALYSIS BY SECTOR (differential image SNR and per-pixel BLS SNR; the full pixel grid is in the data folder)", fontsize=6.5,
                 color=MUTED, fontweight="bold")
        shown = images[:6]
        for i, path in enumerate(shown):
            ax = _axes(fig, .03 + (i % 3) * .325, .21 - (i // 3) * .17, .30, .14)
            image = _trim_white(Image.open(path))
            image = image.crop((0, int(image.height * (1 - PIXEL_MAPS_FRACTION)), image.width, image.height))
            ax.imshow(image)
            ax.axis("off")
            match = re.search(r"folded_tpf_(\d+)", os.path.basename(path))
            ax.set_title(f"sector {match.group(1)}" if match else "", loc="left", fontsize=7, pad=2)
        if len(images) > len(shown):
            fig.text(.97, .385, f"{len(images) - len(shown)} more sectors in the data folder", fontsize=6.5, color=MUTED,
                     ha="right")
    fig.text(.03, .018, "Pink bands: transit duration. Ellipses: 1σ. Hatched: aperture.", fontsize=6.5, color=MUTED)
    return fig


# ---------------------------------------------------------------------------------------------- page 4
def page_transits(inputs, frame, z, flagged):
    """Every measured single transit as a thumbnail; the flagged ones are framed."""
    fig = _new_figure()
    n = len(frame)
    _head(fig, inputs, f"Single transits · {n} measured",
          "Each thumbnail: binned flux (pink: transit duration) and depth in ppt. Framed: the transits with a sheet", 4)
    curve = _lc(inputs)
    mean = np.average(frame.depth, weights=1 / frame.depth_err ** 2)
    a = _axes(fig, .05, .795, .52, .105)
    a.errorbar(frame.index, frame.depth, frame.depth_err, fmt="o", ms=2, color=INK, lw=.4)
    a.axhline(mean, color=OK, lw=.8, ls="--")
    a.set_ylabel("depth (ppt)")
    a.set_xlabel("transit number")
    a.grid(color=GRID, lw=.4)
    _title(a, "Depth of each transit")
    a = _axes(fig, .64, .795, .14, .105)
    a.hist(z, bins=min(14, max(3, n // 2)), color=BLUE, alpha=.8)
    a.set_xlabel("(depth − mean) / error")
    _title(a, "Consistency", f"χ²/n = {(z ** 2).sum() / max(n - 1, 1):.2f}")
    a = _axes(fig, .83, .795, .14, .105)
    a.axis("off")
    listed = len(inputs.transits["t0"]) if inputs.transits is not None else n
    for i, (k, v) in enumerate([("Epochs measured", f"{n} of {listed}"), ("Mean depth", f"{mean:.2f} ppt"),
                                ("Outliers > 2.5σ", f"{int((np.abs(z) > OUTLIER_SIGMA).sum())}"),
                                ("Sheets that follow", f"{min(len(flagged), MAX_SHEETS)}")]):
        a.text(0, .9 - i * .24, k, fontsize=7, color=MUTED, va="top")
        a.text(.62, .9 - i * .24, v, fontsize=8, fontweight="bold", va="top")
    if curve is None:
        _unavailable(_axes(fig, .03, .05, .94, .7), "light curve not available")
        return fig
    t, f = curve
    shown = frame.index[:THUMBNAILS_MAX]
    cols = 12 if len(shown) > 24 else 6
    rows = int(np.ceil(len(shown) / cols))
    x0, y0, w, h = .03, .045, .94, .675
    cw, chh = w / cols, h / max(rows, 4 if cols == 6 else rows)
    dur_d = inputs.duration_h / 24
    for i in shown:
        r = frame.loc[i]
        cx, cy = x0 + (i % cols) * cw, y0 + h - (i // cols + 1) * chh
        ax = _axes(fig, cx + cw * .06, cy + chh * .07, cw * .88, chh * .72)
        lo, hi = np.searchsorted(t, [r.t0 - 3 * dur_d, r.t0 + 3 * dur_d])
        xx = (t[lo:hi] - r.t0) * 24
        if xx.size > 8:
            bx, by, _ = binned(xx, f[lo:hi], xx.min(), xx.max(), 12, minimum=3)
            ax.plot(bx, by, "-", color=INK, lw=.6)
        _band(ax, inputs.duration_h / 2)
        depth = inputs.depth_ppt / 1000
        ax.set_ylim(1 - 2.6 * depth, 1 + 2.6 * depth)
        ax.set_xticks([])
        ax.set_yticks([])
        framed = i in flagged
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color((BAD if abs(z[i]) > OUTLIER_SIGMA else BLUE) if framed else GRID)
            spine.set_linewidth(1.5 if framed else .4)
        fig.text(cx + cw * .08, cy + chh * .83, f"#{i}", fontsize=5.5, color=MUTED)
        fig.text(cx + cw * .95, cy + chh * .83, f"{r.depth:.2f}", fontsize=5.5, ha="right")
    if len(frame) > THUMBNAILS_MAX:
        fig.text(.97, .015, f"First {THUMBNAILS_MAX} of {len(frame)} transits shown", fontsize=6.5, color=MUTED, ha="right")
    fig.text(.03, .015, "Blue frame: deepest, shallowest and closest to the candidate depth. Red frame: |depth − mean| > 2.5 σ.",
             fontsize=6.5, color=MUTED)
    return fig


# ---------------------------------------------------------------------------------------------- page 5+
def _sheet(fig, inputs, tpfs, lc_data_norm, apertures, t, f, row, index, reason, z, y0, hgt):
    """One detailed transit sheet: pixel grid on the left and six panels on the right."""
    from scipy import ndimage
    from watson.watson import Watson
    t0 = float(row.t0)
    tpf = next((p for p in tpfs if p.time.value.min() <= t0 <= p.time.value.max()), None)
    fig.text(.03, y0 + hgt - .003, f"Transit #{index} · T0 {t0:.3f} · depth {row.depth:.2f} ± {row.depth_err:.2f} ppt "
                                  f"({z:+.1f} σ) · {reason}", fontsize=8, fontweight="bold", va="top")
    half_h = 4 * inputs.duration_h
    pw, ph_ = .205, (hgt - .06) / 2 - .012

    def panel(col, r, title):
        ax = _axes(fig, .30 + col * (pw + .015), y0 + .01 + (1 - r) * (ph_ + .03), pw, ph_)
        ax.grid(color=GRID, lw=.3)
        ax.tick_params(labelsize=5.5)
        _title(ax, title)
        return ax

    sel = np.abs(t - t0) < half_h / 24
    xx = (t[sel] - t0) * 24
    depth = inputs.depth_ppt / 1000
    ax = panel(0, 0, "Photometry")
    ax.scatter(xx, f[sel], s=1, color=MUTED, alpha=.35, lw=0)
    if xx.size > 8:
        b = binned(xx, f[sel], -half_h, half_h, 24)
        ax.errorbar(b[0], b[1], b[2], fmt="o", ms=2, color=INK, lw=.5)
    ax.plot(np.linspace(-half_h, half_h, 300), box_model(np.linspace(-half_h, half_h, 300), depth, inputs.duration_h),
            color=BAD, lw=1)
    _band(ax, inputs.duration_h / 2)
    ax.set_ylim(1 - 2.6 * depth, 1 + 2.6 * depth)
    zoom = None
    if lc_data_norm is not None:
        zoom = lc_data_norm[(lc_data_norm["time"] > t0 - half_h / 24) & (lc_data_norm["time"] < t0 + half_h / 24)]
        quality = zoom["quality"].to_numpy() if "quality" in zoom else None
        if quality is not None and not pd.isna(quality).all():
            from watson import constants
            for tm in zoom["time"][np.bitwise_and(np.nan_to_num(quality).astype(int), constants.MOMENTUM_DUMP_QUALITY_FLAG) >= 1]:
                ax.axvline((tm - t0) * 24, color="purple", alpha=.3, ls="--", lw=1)
    for col, r, title, (a_name, a_color, a_label), (b_name, b_color, b_label) in (
            (1, 0, "X-axis drift", ("motion_x", RED, "motion"), ("centroids_x", INK, "centroid")),
            (2, 0, "Y-axis drift", ("motion_y", RED, "motion"), ("centroids_y", INK, "centroid"))):
        ax = panel(col, r, title)
        if zoom is not None and a_name in zoom and len(zoom):
            ax.scatter((zoom["time"] - t0) * 24, zoom[a_name], s=2, color=a_color, label=a_label)
            ax.scatter((zoom["time"] - t0) * 24, zoom[b_name], s=2, color=b_color, label=b_label)
            ax.legend(frameon=False, fontsize=5.5, loc="upper left")
        else:
            _unavailable(ax)
        _band(ax, inputs.duration_h / 2)
    ax_sap, ax_ap, ax_bg = panel(0, 1, "SAP: used vs smaller aperture"), panel(1, 1, "Apertures"), panel(2, 1, "Background flux")
    if zoom is not None and "background_flux" in zoom and len(zoom):
        ax_bg.scatter((zoom["time"] - t0) * 24, zoom["background_flux"], s=2, color=BLUE)
        _band(ax_bg, inputs.duration_h / 2)
    else:
        _unavailable(ax_bg)
    grid_ax = (.03, y0 + .01, .24, hgt - .045)
    if tpf is None or apertures is None:
        _unavailable(ax_sap)
        _unavailable(ax_ap)
        _unavailable(_axes(fig, *grid_ax), "no target pixel file for this transit")
        return
    author, cadence = Watson.get_author_cadence_from_tpf_name(tpf)
    aperture = Watson.get_aperture_for_sector(apertures, _sector_of(tpf), author, cadence)
    from lcbuilder.photometry.aperture_extractor import ApertureExtractor
    mask = ApertureExtractor.from_pixels_to_boolean_mask(aperture, tpf.column, tpf.row, tpf.shape[2], tpf.shape[1])
    eroded = ndimage.binary_erosion(mask)
    window = (tpf.time.value > t0 - half_h / 24) & (tpf.time.value < t0 + half_h / 24)
    flux = tpf.flux.value[window]
    tt = tpf.time.value[window]
    x = (tt - t0) * 24
    median_image = np.nanmedian(tpf.flux.value, axis=0)
    for ap, color, label in ((mask, INK, "used aperture"), (eroded, ORANGE, "smaller aperture")):
        if ap.any() and x.size > 8:
            series = np.nansum(flux[:, ap], axis=1)
            series = series / np.nanmedian(series)
            b = binned(x, series, x.min(), x.max(), 24, True, 2)
            ax_sap.plot(b[0], b[1], "-", color=color, lw=1, label=label)
    _band(ax_sap, inputs.duration_h / 2)
    ax_sap.legend(frameon=False, fontsize=5.5, loc="lower left")
    ny, nx = median_image.shape
    cy, cx = np.unravel_index(np.nanargmax(np.where(mask, median_image, np.nan)), median_image.shape)
    ax_ap.imshow(median_image, origin="lower", cmap="viridis")
    ax_ap.contour(mask, levels=[.5], colors=[RED], linewidths=1)
    if eroded.any():
        ax_ap.contour(eroded, levels=[.5], colors=["#ffd166"], linewidths=1)
    ax_ap.set_xlim(cx - 4.5, cx + 4.5)
    ax_ap.set_ylim(cy - 4.5, cy + 4.5)
    ax_ap.set_xticks([])
    ax_ap.set_yticks([])
    k = 7
    gx0, gy0, gw, gh = grid_ax
    for i in range(k):
        for j in range(k):
            yy, xx_ = cy - k // 2 + i, cx - k // 2 + j
            if not (0 <= yy < ny and 0 <= xx_ < nx):
                continue
            pax = _axes(fig, gx0 + j * gw / k, gy0 + i * gh / k, gw / k * .94, gh / k * .92)
            series = flux[:, yy, xx_]
            med = np.nanmedian(series)
            if x.size > 8 and np.isfinite(med) and med != 0:
                b = binned(x, series / med, x.min(), x.max(), 14, True, 2)
                pax.plot(b[0], b[1], color=INK, lw=.5)
            _band(pax, inputs.duration_h / 2)
            pax.set_xticks([])
            pax.set_yticks([])
            for spine in pax.spines.values():
                spine.set_visible(True)
                spine.set_color(RED if mask[yy, xx_] else GRID)
                spine.set_linewidth(1.3 if mask[yy, xx_] else .3)


def pages_sheets(inputs, frame, z, flagged, first_page):
    """Detailed sheets of the flagged transits, two per page."""
    curve = _lc(inputs)
    if curve is None or not flagged:
        return []
    t, f = curve
    _, _, tpfs, lc_data_norm = _load_tpfs(inputs)
    apertures = _load_apertures(inputs)
    chosen = list(flagged.items())[:MAX_SHEETS]
    pages = []
    for start in range(0, len(chosen), 2):
        fig = _new_figure()
        _head(fig, inputs, "Detailed sheets of the flagged transits · two per page",
              "Photometry, drift, apertures, background and the pixels of every flagged transit", first_page + len(pages))
        for slot, (index, reason) in enumerate(chosen[start:start + 2]):
            try:
                _sheet(fig, inputs, tpfs, lc_data_norm, apertures, t, f, frame.loc[index], index, reason, z[index],
                       .49 if slot == 0 else .03, .43)
            except Exception:
                logging.exception("Could not draw the sheet of transit %s", index)
        if start == 0 and len(flagged) > MAX_SHEETS:
            fig.text(.97, .012, f"{len(flagged) - MAX_SHEETS} more flagged transits without a sheet", fontsize=6.5,
                     color=MUTED, ha="right")
        pages.append(fig)
    return pages


# ---------------------------------------------------------------------------------------------- writer
def _stub(inputs, title, error, page):
    fig = _new_figure()
    _head(fig, inputs, title, "This page could not be built", page)
    fig.text(.5, .5, f"Not available: {error}", ha="center", color=MUTED, fontsize=9)
    return fig


def build_pages(inputs):
    """All the pages of the report; a page that fails is replaced by a stub that says why."""
    metrics = load_metrics(inputs.data_dir)
    pages = []

    def add(title, builder):
        try:
            with matplotlib.rc_context(RC):
                built = builder(len(pages) + 1)
        except Exception as error:
            logging.exception("The page '%s' of the vetting report failed", title)
            built = [_stub(inputs, title, f"{type(error).__name__}: {error}", len(pages) + 1)]
        pages.extend(built if isinstance(built, list) else [built])

    add("Verdict and signal", lambda n: page_verdict(inputs, metrics))
    add("Classifiers", lambda n: page_classifiers(inputs, metrics))
    add("Is the signal on the target?", lambda n: page_on_target(inputs, metrics))
    if inputs.transits is not None:
        frame, z, flagged = transit_depths_summary(inputs.transits, inputs.depth_ppt)
        if len(frame):
            add("Single transits", lambda n: page_transits(inputs, frame, z, flagged))
            add("Flagged transits", lambda n: pages_sheets(inputs, frame, z, flagged, n))
    return pages


def write_landscape_report(inputs, path):
    """Write the report atomically (``.partial`` first) and return the number of pages."""
    pages = build_pages(inputs)
    partial = path + ".partial"
    with matplotlib.rc_context({"pdf.fonttype": 42}), PdfPages(partial) as pdf:  # TrueType: searchable text
        for fig in pages:
            pdf.savefig(fig)
    os.replace(partial, path)
    return len(pages)
