#!/bin/bash

rm tests.log
rm dist* -r
set -e
tox -r > tests.log
tests_results=$(cat tests.log | grep "congratulations")
if ! [[ -z ${tests_results} ]]; then
  git_tag=$1
  sed -i '6s/.*/version = "'${git_tag}'"/' setup.py
  sed -i '1s/.*/__version__ = "'${git_tag}'"/' watson/__init__.py
  git add setup.py
  git add watson/__init__.py
  git commit -m "Preparing release ${git_tag}"
  git tag ${git_tag} -m "New release"
  git push
  git push --tags
  python3 setup.py sdist bdist_wheel
  python3 -m twine upload dist/*
  rm tests.log
else
  echo "TESTS FAILED. See tests.log"
fi
rm dist* -r
rm -r .tox

