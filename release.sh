#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate base
rm tests.log
rm dist* -r
rm -r .pytest_cache
rm -r build
rm -r watson-reqs
conda remove -n watson-reqs --all -y
rm -R *egg-info
rm -R .tox
set -e
tox > tests.log
tests_results=$(cat tests.log | grep "congratulations")
if ! [[ -z ${tests_results} ]]; then
  echo "Building"
  set +e
  rm tests.log
  rm dist* -r
  rm -r .pytest_cache
  rm -r build
  rm -r watson-reqs
  conda remove -n watson-reqs --all -y
  set -e
  conda create -n watson-reqs python=3.13 -y
  conda activate watson-reqs
  python3 -m pip install pip -U
  python3 -m pip install numpy==2.2.4
  git_tag=$1
  sed -i '6s/.*/version = "'${git_tag}'"/' setup.py
  sed -i '1s/.*/__version__ = "'${git_tag}'"/' watson/__init__.py
  python3 -m pip install .
  python3 -m pip list --format=freeze > requirements.txt
  conda deactivate
  git pull
  git add setup.py
  git add watson/__init__.py
  git add requirements.txt
  git commit -m "Preparing release ${git_tag}"
  git tag ${git_tag} -m "New release"
  git push && git push --tags
else
  echo "Failed tests"
fi
set +e
rm -R watson-reqs
rm dist* -r
rm -r .pytest_cache
rm -r build
rm -R *egg-info
rm -R .tox
conda remove -n watson-reqs --all -y
set -e
