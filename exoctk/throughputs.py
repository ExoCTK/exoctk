# !/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for creating and managing grids of model spectra
"""

from glob import glob
import json
import os
from pathlib import Path

import numpy as np
from svo_filters.svo import Filter

from .pkgdata import resource_filename

try:
    from pandeia.engine.calc_utils import build_default_calc
    from pandeia.engine.custom_exceptions import EngineInputError
    from pandeia.engine.instrument_factory import InstrumentFactory
except ImportError:
    print("pandeia not installed. Functionality limited.")

def get_throughput_list(filt_dir):
    filt_dir = Path(filt_dir)
    filt_files = [x.stem for x in filt_dir.rglob('*.txt')]
    return filt_files

FILT_DIR = resource_filename('exoctk', 'data/throughputs/')
JWST_THROUGHPUTS = get_throughput_list(FILT_DIR)


class Throughput(Filter):
    def __init__(self, name, **kwargs):
        """
        Initialize the Throughput object as a child class of
        svo_filters.svo.Filter

        Parameters
        ----------
        name: str
            The [instrument.element_name] of the filter or disperser
        """
        # Check if JWST instrument
        if name in JWST_THROUGHPUTS:

            # Change throughput directory
            super().__init__(name, filter_directory=Path(FILT_DIR), **kwargs)

        # Otherwise, just use svo_filters default
        else:
            super().__init__(name, **kwargs)


def get_pce(instrument='niriss', mode='soss', filter='clear', disperser='gr700xd', aperture='soss', detector=None, order=None):
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
    detector: dict (optional)
        The detector configuration to use. Defaults to the pandeia default
        detector configuration for the given instrument and mode.
    order: int (optional)
        The spectral order to use for multi-order dispersers (e.g. NIRISS SOSS)

    Returns
    -------
    wave, pce
        The wavelength and throughput
    """
    # Start from the pandeia default configuration so mode-specific settings
    # (e.g. a valid detector configuration) are populated
    conf = build_default_calc('jwst', instrument, mode)['configuration']
    conf['instrument'].update({'filter': filter, 'aperture': aperture, 'disperser': disperser})
    if detector is not None:
        conf['detector'] = detector

    # NIRCam coronagraphic apertures encode the detector channel (e.g. 'mask335rlw')
    if 'detector' in conf['instrument'] and aperture[-2:] in ('sw', 'lw'):
        conf['instrument']['detector'] = aperture[-2:]

    i = InstrumentFactory(config=conf)
    i.order = order

    # Determine wavelength range
    wr = i.get_wave_range()
    wave = np.linspace(wr['wmin'], wr['wmax'], num=500)

    # Evaluate the throughput, clipping small negative interpolation artifacts
    pce = np.clip(i.get_total_eff(wave), 0, None)

    return wave, pce


def _skip_combination(instrument, mode, aperture, filt, disperser=None):
    """
    Determine whether a pandeia configuration should be excluded from the
    generated throughput files

    Parameters
    ----------
    instrument: str
        The JWST instrument
    mode: str
        The observing mode
    aperture: str
        The aperture
    filt: str
        The filter
    disperser: str
        The disperser

    Returns
    -------
    bool
        True if the combination should be skipped
    """
    # Only keep NIRCam coronagraphic masks on their native detector channel
    if instrument == 'nircam' and mode == 'coronagraphy':
        native = 'sw' if aperture[:-2] in ('mask210r', 'maskswb') else 'lw'
        if aperture[-2:] != native:
            return True

        # Pandeia allows F277W with the round LW occulters, but JDox does not
        # (F277W is bar-only and F322W2 is round-only)
        return filt == 'f277w' and aperture[:-2] in ('mask335r', 'mask430r')

    # GRISMC has the same throughput as GRISMR and is not used for time series
    if instrument == 'nircam' and disperser == 'grismc':
        return True

    # MIRI WFSS uses the same P750L + imager configuration as LRS slitless
    if instrument == 'miri' and mode == 'wfss':
        return True

    # MIRI coronagraphic filters are only used with their coronagraphic masks
    if instrument == 'miri' and mode != 'coronagraphy' and filt is not None and filt.endswith('c'):
        return True

    return False


def generate_JWST_throughputs(path=None, data_dir=None):
    """
    Function to generate .txt files of all JWST filter and grism throughputs

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

        # Iterate over science instruments (NIRSpec separated out below)
        for inst in ['niriss', 'nircam', 'miri']:

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
                    filt_list = params['filters'] or [None]
                    disp_list = params.get('dispersers') or [None]

                    # Make a list of all aperture, disperser, order, and filter combinations,
                    # splitting dispersers into orders if necessary
                    combinations = []
                    for aper in aper_list:
                        for disp in disp_list:
                            for filt in filt_list:
                                orders = cfg['disperser_config'].get(disp, {}).get('orders', [1]) if disp else [1]
                                if isinstance(orders, dict):
                                    orders = orders.get(filt, [1])
                                for order in (orders if len(orders) > 1 else [None]):
                                    if not _skip_combination(inst, mode, aper, filt, disp):
                                        combinations.append([aper, disp, order, filt])

                    for conf in combinations:

                        try:

                            # Get the configuration
                            aper, disp, order, filt = conf

                            # Generate throughput
                            wave, thru = get_pce(inst, mode, filt, disp, aper, order=order)
                            disp_name = 'CLEARP' if disp is None else disp if order is None else '{}_{}'.format(disp, order)
                            aper_name = aper[:-2] if inst == 'nircam' and mode == 'coronagraphy' else aper
                            if inst == 'miri' and mode == 'lrsslitless':
                                aper_name = 'lrsslitless'
                            name = '.'.join([inst.upper().replace('CAM', 'Cam').replace('SPEC', 'Spec'), (filt or 'CLEAR').upper(), disp_name.upper(), aper_name.upper()])

                            # Save to txt file
                            data = np.array([wave, thru]).T
                            datafile_path = os.path.join(data_dir, '{}.txt'.format(name))
                            with open(datafile_path, 'w+') as datafile_id:
                                np.savetxt(datafile_id, data)
                            print("{} file created!".format(datafile_path))

                        except (KeyError, ValueError, EngineInputError) as err:
                            print("Could not produce throughput for {}: {}".format(conf, err))

        # Do NIRSpec separately because pandeia doesn't show them correctly
        nirspec_combos = ['G140M/F070LP', 'G140M/F100LP', 'G235M/F170LP', 'G395M/F290LP', 'G140H/F070LP', 'G140H/F100LP', 'G235H/F170LP', 'G395H/F290LP', 'PRISM/CLEAR']

        for disp, filt in [i.split('/') for i in nirspec_combos]:
            inst = 'nirspec'
            aper = 's200a1'
            mode = 'fixed_slit'
            detector = {"nexp": 1, "ngroup": 10, "nint": 1, "readout_pattern": "nrs", "subarray": "full"}

            # Generate throughput
            wave, thru = get_pce(inst, mode, filt.lower(), disp.lower(), aper, detector)
            name = '.'.join(
                [inst.upper().replace('CAM', 'Cam').replace('SPEC', 'Spec'), filt.upper(),
                 disp.upper(), aper.upper()])

            # Save to txt file
            data = np.array([wave, thru]).T
            datafile_path = os.path.join(data_dir, '{}.txt'.format(name))
            with open(datafile_path, 'w+') as datafile_id:
                np.savetxt(datafile_id, data)
            print("{} file created!".format(datafile_path))
