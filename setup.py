#!/usr/bin/env python
from setuptools import setup

setup(name='ExoCTK',
      version=0.2,
      description='Stuff',
      install_requires=['numpy', 'astropy', 'scipy', 'cython', 'matplotlib', 'numba', 'pysynphot', 'sphinx_automodapi', 'sphinx_rtd_theme', 'bibtexparser', 'bokeh', 'batman-package', 'pandas', 'lmfit', 'svo_filters', 'sphinx_astropy'],
      author='Joe Filippazzo and Jonathan Fraine',
      author_email='jfilippazzo@stsci.edu',
      license='MIT',
      url='https://github.com/spacetelescope/awesimsoss',
      long_description='More stuff',
      zip_safe=False,
      use_2to3=False
)
