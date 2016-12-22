#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module of plotting tools for the limb darkening subpackage.
"""
import matplotlib
import matplotlib.pyplot as plt

def ld_v_mu(coeffs, mu):
    """
    Plot the limb darkening versus mu
    
    Parameters
    ----------
    
    """
    # Create the plot
    plt.figure()
    
    # 