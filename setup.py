#!/usr/bin/env python
import os
from setuptools import setup, find_packages

REQUIRES = ['asteval',
            'astropy',
            'astroquery',
            'bandit',
            'bibtexparser',
            'bokeh',
            'boto3',
            'corner',
            'cython',
            'docopt',
            'docutils',
            'flake8',
            'flask',
            'flask_wtf',
            'gunicorn',
            'h5py',
            'hotsoss',
            'ipython',
            'matplotlib',
            'numpy',
            'numpydoc',
            'pandas',
            'paramiko',
            'platon',
            'pysiaf',
            'pysynphot',
            'pytest',
            'pyyaml',
            'pyvo',
            'regions',
            'scipy',
            'scp',
            'sphinx',
            'sphinx_astropy',
            'sqlalchemy',
            'svo_filters',
            'wtforms',
            'werkzeug',
            'jwst_gtvt']


FILES = []
for root, _, files in os.walk("exoctk"):
    FILES += [os.path.join(root.replace("exoctk/", ""), fname)
              for fname in files if not fname.endswith(".py") and not fname.endswith(".pyc")]

setup(
    name='exoctk',
    version='1.2.5',
    description='Observation reduction and planning tools for exoplanet science',
    packages=find_packages(
        ".",
        exclude=["*.tests"]),
    package_data={
        'exoctk': FILES},
    install_requires=REQUIRES,
    author='The ExoCTK Group',
    author_email='exoctk@gmail.com',
    license='MIT',
    url='https://github.com/ExoCTK/exoctk',
    long_description='',
    zip_safe=True,
    use_2to3=False)
