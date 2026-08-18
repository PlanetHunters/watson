#!/bin/bash

set -e

run_tests=false
git_tag=""
for arg in "$@"; do
  if [[ "${arg}" == "--tests" ]]; then
    run_tests=true
  elif [[ -z "${git_tag}" ]]; then
    git_tag="${arg}"
  else
    echo "Unknown argument: ${arg}" >&2
    exit 1
  fi
done
if [[ -z "${git_tag}" ]]; then
  echo "Usage: $0 VERSION [--tests]" >&2
  exit 1
fi

CONDA_EXE="${CONDA_EXE:-$(command -v conda || true)}"
if [[ -z "${CONDA_EXE}" ]]; then
  echo "Could not find conda. Set CONDA_EXE before running this script." >&2
  exit 1
fi
CONDA_BASE="$(${CONDA_EXE} info --base)"
export CONDA_EXE
export PATH="${CONDA_BASE}/condabin:${CONDA_BASE}/bin:${PATH}"

rm -f tests.log
rm -rf dist* .tox .pytest_cache build watson-reqs *egg-info
"${CONDA_EXE}" env remove -n watson-reqs -y 2>/dev/null || true
if [[ "${run_tests}" == true ]]; then
  if python -m tox -r > tests.log; then
    tests_results=$(grep "congratulations" tests.log || true)
    if [[ -z "${tests_results}" ]]; then
      echo "Tests failed; release aborted." >&2
      exit 1
    fi
  else
    echo "Tests failed; release aborted." >&2
    exit 1
  fi
else
  echo "Skipping tests (use --tests to run them)."
fi

"${CONDA_EXE}" create --override-channels -c conda-forge -n watson-reqs python=3.11 -y
RELEASE_PYTHON="${CONDA_BASE}/envs/watson-reqs/bin/python"
"${RELEASE_PYTHON}" -m pip install -U pip setuptools numpy==2.1.1
sed -i '6s/.*/version = "'${git_tag}'"/' setup.py
sed -i '1s/.*/__version__ = "'${git_tag}'"/' watson/__init__.py
"${RELEASE_PYTHON}" -m pip install .
"${RELEASE_PYTHON}" -m pip list --format=freeze > requirements.txt
git pull
git add setup.py watson/__init__.py requirements.txt
git commit -m "Preparing release ${git_tag}"
git tag "${git_tag}" -m "New release"
git push
git push --tags
set +e
rm -rf watson-reqs dist* .tox .pytest_cache build *egg-info
"${CONDA_EXE}" env remove -n watson-reqs -y 2>/dev/null || true
