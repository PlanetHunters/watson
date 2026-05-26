import logging
import os
import sys


def _setup_logging():
    log_dir = os.environ.get("PYTEST_LOG_DIR", os.path.dirname(__file__))
    log_file = os.path.join(log_dir, "pytest_output.log")

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)-8s] %(name)s %(filename)s:%(lineno)d %(message)s",
        datefmt="%H:%M:%S",
    )

    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(fmt)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.addHandler(fh)
    root.addHandler(sh)


_setup_logging()


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
