#!/usr/bin/env python
import os
from setuptools import setup, find_packages

REQUIRES = ['asteval',
            'astropy',
            'astroquery',
            'bibtexparser',
            'bokeh',
            'cython',
            'flask',
            'h5py',
            'lmfit',
            'matplotlib',
            'numba',
            'numpy',
            'pandas',
            'pysynphot',
            'scipy',
            'sphinx',
            'svo_filters']

FILES = []
for root, _, files in os.walk("exoctk"):
    FILES += [os.path.join(root.replace("exoctk/", ""), fname) \
        for fname in files if not fname.endswith(".py") and not fname.endswith(".pyc")]

setup(name='exoctk',
      version='0.2.2',
      description='Observation reduction and planning tools for exoplanet science',
      packages=find_packages(".", exclude=["*.tests"]),
      package_data={'exoctk': FILES},
      install_requires=REQUIRES,
      author='The ExoCTK Group',
      author_email='exoctk@gmail.com',
      license='MIT',
      url='https://github.com/ExoCTK/exoctk',
      long_description='',
      zip_safe=True,
      use_2to3=False
)
