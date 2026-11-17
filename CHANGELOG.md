# Changelog

## Unreleased

### Changed

- The vetting report is now a compact landscape PDF: verdict and signal,
  WATSON-NET and TRICERATOPS, source offset / centroids / optical ghost with
  the per-sector pixel maps, every single transit as a thumbnail, and detailed
  sheets only for the flagged transits (deepest, shallowest, closest to the
  candidate depth and every depth outlier above 2.5 sigma, at most 12).
  A TIC with 62 transits goes from 83 pages / 66 MB to 7 pages / under 1 MB.
- The separate `*_summary.pdf` is no longer generated and `only_summary` is
  ignored by the landscape report. The per-transit PNGs are only produced with
  `legacy_report=True`.
- If the landscape report cannot be built, `Watson.vetting` falls back to the
  legacy portrait report. A page that fails becomes a stub that says why.

### Added

- `legacy_report` argument of `Watson.vetting` and `Watson.vetting_with_data`
  to keep the old portrait report.

## 2.0.0 - 2026-10-17

### Breaking changes

- WATSON-NET scores now weight TRICERATOPS planet and false-positive scenarios
  directly; the previous training-set prior correction is removed.
- `Watson.compute_bayesian_fpp` no longer accepts
  `training_planet_fraction`, and `bayesian_fpp.csv` no longer includes that
  column.

### Added

- Forward optional external ground-based light curves and filter names to
  TRICERATOPS+.
- Add an integration-test group covering network-backed vetting and report
  generation.

### Fixed and updated

- Pin `MeepMeep==0.8.0` to remain compatible with `PyTransit==2.6.14`.
- Update the ExoML dependency to 1.4.0.
- Use the selected TRICERATOPS sigma mode directly and remove the legacy
  Eleanor fallback from neighbor-lightcurve processing.
- Improve MAST TLS handling, report generation, and API documentation.
