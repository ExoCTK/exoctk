# !/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for replacing the no-longer-available pkg_resources module in python 3.12+
"""

from importlib import resources


def resource_filename(package_name, resource_name):
    """Name-alike replacement for pkg_resources.resource_filename"""
    # Use a context manager to ensure the resource path is valid, 
    # especially if the package is in a zip archive.
    ref = resources.files(package_name).joinpath(resource_name)
    with resources.as_file(ref) as path:
        # path is a pathlib.Path object pointing to the file
        # You can now use this path as a string (e.g., in open())
        return str(path)
