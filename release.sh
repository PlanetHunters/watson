#!/bin/bash

git_tag=$1
sed -i '6s/.*/version = "'${git_tag}'"/' setup.py
sed -i '1s/.*/__version__ = "'${git_tag}'"/' watson/__init__.py
git add setup.py
git add watson/__init__.py
git commit -m "[ci skip] Preparing release ${git_tag}"
git tag ${git_tag} -m "New release"
git push -o ci.skip && git push --tags
