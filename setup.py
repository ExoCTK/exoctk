#!/usr/bin/env python
import os
from setuptools import setup, find_packages

REQUIRES = ['astropy==5.1',
            'bokeh==2.4.3',
            'boto3',
            'cython==0.29.32',
            'docopt==0.6.2',
            'docutils==0.16.0',
            'flake8==4.0.1',
            'flask==2.1.3',
            'gunicorn==20.1.0',
            'h5py==3.7.0',
            'ipython==8.6.0',
            'jinja2==3.1.2',
            'jupyter==1.0.0',
            'matplotlib==3.6.2',
            'numpy==1.23.4',
            'numpydoc==1.5.0',
            'pandas==1.5.2',
            'paramiko==2.8.1',
            'pip==22.2.2',
            'pytest==7.1.2', 
            'pyyaml==5.4',
            'scipy==1.9.3',
            'scp==0.14.1',
            'sphinx==5.0.2',
            'sqlalchemy==1.4.39',
            'twine==3.7.1',
            'wtforms==2.3.3',
            'asteval==0.9.28',
            'awscli',
            'bandit==1.7.4',
            'batman-package==2.4.9',
            'bibtexparser==1.4.0',
            'corner==2.2.1',
            'ddtrace==1.6.3',
            'flask_wtf==1.0.1',
            'lmfit==1.1.0',
            'platon==5.4',
            'pysiaf==0.19.0',
            'pysynphot==2.0.0',
            'sphinx-astropy==1.7.0',
            'svo-filters==0.4.4',
            'werkzeug==2.1.1',
            'jwst_gtvt==0.3.1',
            'astroquery==0.4.6']

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
