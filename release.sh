#!/bin/bash


rm tests.log
rm dist* -r
rm -r .tox
rm -r .pytest_cache
rm -r build
set -e
tox -r > tests.log
tests_results=$(cat tests.log | grep "congratulations")
if ! [[ -z ${tests_results} ]]; then
  python3.8 -m venv watson-reqs
  source watson-reqs/bin/activate
  python3.8 -m pip install pip -U
  python3.8 -m pip install numpy==1.22.3
  git_tag=$1
  sed -i '6s/.*/version = "'${git_tag}'"/' setup.py
  sed -i '1s/.*/__version__ = "'${git_tag}'"/' watson/__init__.py
  python3.8 setup.py install
  python3.8 -m pip list --format=freeze > requirements.txt
  deactivate
  git pull
  git add setup.py
  git add watson/__init__.py
  git add requirements.txt
  git commit -m "Preparing release ${git_tag}"
  git tag ${git_tag} -m "New release"
  git push && git push --tags
fi
rm -R watson-reqs
rm dist* -r
rm -r .tox
rm -r .pytest_cache
rm -r build