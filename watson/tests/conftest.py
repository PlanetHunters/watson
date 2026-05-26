# conftest.py for watson tests

import logging
import sys


def _setup_logging():
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)-8s] %(filename)s:%(lineno)d %(message)s",
        datefmt="%H:%M:%S",
    ))
    root.handlers.clear()
    root.addHandler(handler)


_setup_logging()


# Ensure pkg_resources is importable before any test modules load.
# pytransit (a transitive dependency via triceratops) imports pkg_resources
# at module load time, and modern setuptools (>=71) no longer provides it.
# This early import forces a clear error if the environment is misconfigured.
try:
    import pkg_resources  # noqa: F401
except ImportError:
    import warnings
    warnings.warn(
        "pkg_resources is not available. Make sure 'setuptools<71' is installed "
        "in the test environment (tox.ini / CI).",
        RuntimeWarning,
    )
    raise
