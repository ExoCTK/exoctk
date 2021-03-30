# !/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for creating and managing grids of model spectra
"""
from glob import glob
import json
import os
from pkg_resources import resource_filename
import warnings

import numpy as np
from svo_filters.svo import Filter

try:
    from pandeia.engine.instrument_factory import InstrumentFactory
except ImportError:
    print("pandeia not installed. Functionality limited.")

FILT_DIR = resource_filename('exoctk', 'data/throughputs/')
JWST_THROUGHPUTS = [os.path.basename(file).replace('.txt', '') for file in glob(FILT_DIR + '*')]


class Throughput(Filter):
    def __init__(self, name, **kwargs):
        """
        Initialize the Throughput object as a child class of svo_filters.svo.Filter

        Parameters
        ----------
        name: str
            The [instrument.element_name] of the filter or disperser
        """
        # Check if JWST instrument
        if name in JWST_THROUGHPUTS:

            # Change throughput directory
            super().__init__(name, filter_directory=FILT_DIR, **kwargs)

        # Otherwise, just use svo_filters default
        else:
            super().__init__(name, **kwargs)


def get_pce(instrument='niriss', mode='soss', filter='clear', disperser='gr700xd', aperture='soss'):
    """
    Use pandeia to generate the JWST throughputs

    Parameters
    ----------
    instrument: str
        The JWST instrument to use, ['niriss', 'nirspec', 'miri', 'nircam']
    mode: str
        The observing mode to use, ['soss', 'wfss', 'ami', etc.]
    filter: str
        The filter wheel element to use
    disperser: str
        The dispersion element to use
    aperture: str
        The aperture to use

    Returns
    -------
    wave, pce
        The wavelength and throughput
    """
    # Get optical element configuration
    obsmode = {'instrument': instrument, 'mode': mode, 'filter': filter, 'aperture': aperture, 'disperser': disperser}
    conf = {'instrument': obsmode}
    i = InstrumentFactory(config=conf)

    # Determine wavelength range
    wr = i.get_wave_range()
    wave = np.linspace(wr['wmin'], wr['wmax'], num=500)

    # Evaluate the throughput
    pce = i.get_total_eff(wave)

    return wave, pce


def generate_JWST_throughputs(path=None, data_dir=None):
    """
    Function to generate .txt filte of all JWST filter and grism throughputs

    Parameters
    ----------
    path: str
        The path to the directory
    data_dir
    """
    # Check if environment variable exists
    path = path or os.environ.get('pandeia_refdata')
    data_dir = data_dir or resource_filename('exoctk', 'data/throughputs/')

    if path is None:
        raise IOError("No path to pandeia_refdata directory provided!")

    else:

        # Iterate over science instruments
        for inst in ['niriss', 'nircam', 'nirspec', 'miri']:

            # Get all valid configurations from pandeia config file
            cfgfile = open(os.path.join(path, 'jwst/{}/config.json'.format(inst)))
            cfg = json.load(cfgfile)
            cfgfile.close()

            # Iterate over list of modes
            for mode, params in cfg['mode_config'].items():

                # Only interested in science modes
                if mode != 'target_acq':

                    # Get aperture and filter lists
                    aper_list = params['apertures']
                    filt_list = params['filters'] or ['CLEAR']

                    # Split disperser into orders if necessary
                    disp_list = []
                    try:
                        for dis in params['dispersers']:
                            orders = cfg['disperser_config'][dis].get('orders', [1])
                            if len(orders) == 1:
                                disp_list.append(dis)
                            else:
                                for order in orders:
                                    disp_list.append('{}_{}'.format(dis, order))
                    except KeyError:
                        disp_list = None
                    disp_list = disp_list or ['CLEARP']

                    # Make a list of all aperture, dispersion, and filter combinations
                    combinations = [[[[a, d, f] for f in filt_list] for d in disp_list] for a in aper_list]
                    combinations = [k for l in [i for j in combinations for i in j] for k in l]

                    for conf in combinations:

                        try:

                            # Get the configuration
                            aper, disp, filt = conf

                            # Generate throughput
                            wave, thru = get_pce(inst, mode, filt, disp, aper)
                            name = '.'.join([inst.upper().replace('CAM', 'Cam').replace('SPEC', 'Spec'), filt.upper(), disp.upper(), aper.upper()])

                            # Save to txt file
                            data = np.array([wave, thru]).T
                            datafile_path = os.path.join(data_dir, '{}.txt'.format(name))
                            with open(datafile_path, 'w+') as datafile_id:
                                np.savetxt(datafile_id, data)
                            print("{} file created!".format(datafile_path))

                        except KeyError:
                            print("Could not produce throughput for {}".format(conf))