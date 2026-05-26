# AGENTS.md — WATSON (dearwatson)

> Agent-facing guide for the `watson` (dearwatson) project.

## Project Intent

**WATSON** (Visual Vetting and Analysis of Transits from Space ObservatioNs) is a Python package for visually vetting transiting exoplanet candidates detected by space-based photometry missions (TESS, Kepler, K2). It automates many checks found in NASA SPOC Data Validation reports, including transit shape analysis, odd-even checks, weak secondary detection, centroid shifts, optical ghost diagnostics, FPP computation via TRICERATOPS, and optional AI analysis (WATSON-NET, GPT-4o).

- **Package:** `dearwatson`
- **Author:** M. Dévora-Pajares
- **Repository:** https://github.com/PlanetHunters/watson
- **Documentation:** https://dearwatson.readthedocs.io

## Architecture

```
watson/
├── watson.py                    # Core class (~2700 lines)
├── report.py                    # PDF report generation
├── neighbours.py                # Neighbour star catalog queries
├── data_validation_report/
│   └── DvrPreparer.py          # Download official DV reports
├── tpfplotterSub/               # TPF/FOV plotting subpackage
│   ├── tpfplotter.py
│   └── tpfplotter_py2.py
├── tests/
│   └── test_watson.py
└── resources/
```

### Key Classes

| Class | Role |
|-------|------|
| `Watson` | Main orchestrator. `vetting(...)` runs the full pipeline. |
| `TriceratopsThreadValidator` | Wraps TRICERATOPS for false-positive validation. |
| `Report` | Generates PDF vetting reports via ReportLab. |
| `DvrPreparer` | Downloads official mission Data Validation Reports. |

### Execution Flow

1. Instantiate `Watson(object_dir, output_dir)`.
2. Call `vetting(id, period, t0, duration, depth, ...)`.
3. Optionally downloads light curves / TPFs via `lcbuilder`.
4. Runs TRICERATOPS validation to compute FPP/NFPP.
5. Generates folded-curve, odd/even, secondary-event plots.
6. Computes per-cadence SNRs, bootstrap FAP.
7. If TPFs exist: runs centroid shifts, optical ghost, source offset, per-pixel BLS in parallel.
8. Optionally runs WATSON-NET and/or GPT-4o analysis.
9. Writes `metrics.csv` and generates PDF report(s).

## Entrypoints

No registered `console_scripts`. Primary usage is as a Python library.

### Library Entrypoint

```python
from watson.watson import Watson
import os

watson = Watson(object_dir='./data', output_dir='./results')
watson.vetting(
    id="TIC 307210830",
    period=1.049,
    t0=1354.716,
    duration=54,       # minutes
    depth=0.158,       # ppt
    depth_err=0.01,
    sectors='all',
    rp_rstar=0.01205,
    a_rstar=9,
    ra=124.531756,
    dec=-68.312999,
    cadence=[120],
    cpus=os.cpu_count() // 2,
    clean=True
)
```
## Environment Setup

### Prerequisites

- Python `>= 3.11`
- `conda` (tox-conda is used for tests)
- OS build tools: `build-essential`, `libssl-dev`, `gfortran`, `gcc`

### Conda Environment Setup

```bash
cd /path/to/watson
conda create -n watson311 python=3.11
conda activate watson311

# Install build deps first (critical pin)
pip install numpy==2.1.1 Cython==3.0.6 wheel 'setuptools<71'

# Then install the package
pip install -e .
```

**Important:** `setuptools` must be `<71` because `pytransit==2.6.14` requires `pkg_resources` at load time.

## Build & Test

### Run Tests (tox — recommended)

```bash
tox
```

- Uses `tox-conda` to provision a Python 3.11 environment.
- Installs `pytest`, `pytest-forked`, `pytest-xdist`, `numpy==2.1.1`, `Cython`, `setuptools<71`.
- Command: `pytest --forked -n 1 -v -x watson/tests/`

### Run Tests (pytest directly)

```bash
pytest --forked -n 1 -v -x watson/tests/
```

## Execution Recipes

### Use Case A: Import and run as a library

```python
import os
from watson.watson import Watson

object_dir = os.getcwd() + '/my_data'
output_dir = os.getcwd() + '/my_results'
os.makedirs(output_dir, exist_ok=True)

watson = Watson(object_dir=object_dir, output_dir=output_dir)
watson.vetting(
    id="TIC 307210830",
    period=1.0491800670761966,
    t0=1354.7155898963902,
    duration=54,
    depth=0.158,
    depth_err=0.01,
    sectors='all',
    rp_rstar=0.01205,
    a_rstar=9,
    ra=124.531756290083,
    dec=-68.3129998725044,
    cadence=120,
    cpus=os.cpu_count() // 2,
    clean=True
)
```

**Result:** `output_dir` will contain PNG plots, CSV metrics, TRICERATOPS results, and PDF reports (`*_transits_validation_report.pdf`, `_summary.pdf`).

### Use Case B: Run the example notebook

```bash
cd /path/to/watson/examples/TOI-175
jupyter notebook TOI-175-vetting.ipynb
```

### Use Case C: Run tests

```bash
cd /path/to/watson
tox
# or directly:
pytest --forked -n 1 -v -x watson/tests/
```
