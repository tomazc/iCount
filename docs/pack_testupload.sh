#!/usr/bin/env bash
pip install twine
pip install pypandoc
python setup.py clean -a
rm dist/*
rm -r *.egg-info
python setup.py sdist
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
