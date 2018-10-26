#!/usr/bin/env python
from setuptools import setup

setup(name='exoctk',
      version=0.2,
      description='Observation reduction and planning tools for exoplanet science',
      install_requires=['numpy', 'astropy', 'scipy', 'cython', 'matplotlib', 'numba', 'pysynphot', 'sphinx_automodapi', 'sphinx_rtd_theme', 'bibtexparser', 'bokeh', 'pandas', 'svo_filters', 'sphinx_astropy', 'batman-package', 'lmfit', 'flask'],
      author='The ExoCTK Group',
      author_email='jfilippazzo@stsci.edu',
      license='MIT',
      url='https://github.com/ExoCTK/exoctk',
      long_description='',
      zip_safe=False,
      use_2to3=False
)
