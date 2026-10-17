# Changelog

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
