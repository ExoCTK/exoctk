#!/usr/bin/env python
import os
from setuptools import setup, find_packages

REQUIRES = ['asteval',
            'astropy<4.1',
            'astroquery',
            'bandit',
            'bibtexparser',
            'bokeh',
            'boto3',
            'corner',
            'cython',
            'docopt',
            'docutils==0.15.2',
            'flake8',
            'flask',
            'flask_wtf',
            'gunicorn',
            'h5py',
            'ipython',
            'lmfit',
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
            'scipy',
            'scp',
            'sphinx',
            'sphinx_astropy',
            'sqlalchemy',
            'svo_filters',
            'wtforms',
            'werkzeug==0.16.1',
            'jwst_gtvt @ git+https://github.com/spacetelescope/jwst_gtvt.git@cd6bc76f66f478eafbcc71834d3e735c73e03ed5']

FILES = []
for root, _, files in os.walk("exoctk"):
    FILES += [os.path.join(root.replace("exoctk/", ""), fname)
              for fname in files if not fname.endswith(".py") and not fname.endswith(".pyc")]

setup(
    name='exoctk',
    version='1.1.1',
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
