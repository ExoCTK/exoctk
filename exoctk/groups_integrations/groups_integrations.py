"""This is a module for calcuating groups and integrations with JWST on
the fly. The main function (``perform_calculation``) takes a dictionary
of inputs (modeled around how the web tool takes inputs) which must
include: ``observation_time``, ``num_groups``, ``magnitude``,
``model``, ``band``, ``filt``, ``filt_ta``, ``instrument``,
``subarray``, ``subarray_ta``, ``saturation_mode``, ``max_saturation``,
and ``infile``. It produces and dictionary of outputs that includes all
of the original information as well as groups, integrations, saturation
levels, and observation time estimates for target acquisition and
science observations with JWST.

Authors
-------

    Jules Fowler, April 2017
    Matthew Bourque, February 2021

Use
---

    This is mostly a module to be used by the ExoCTK web application,
    but the main function can be run standalone as such:
    ::

        from exoctk.groups_integrations.groups_integrations import perform_calculation
        perform_calculation()

Dependencies
------------

    - ``astropy``
    - ``numpy``
    - ``scipy``
"""

import json
import math

import numpy as np
from scipy import interpolate


def calc_duration_time(num_groups, num_integrations, num_reset_frames, frame_time, frames_per_group=1):
    """Calculates duration time (or exposure duration as told by APT)

    Parameters
    ----------
    num_groups : int
        Groups per integration
    num_integrations : int
        Integrations per exposure
    num_reset_frames : int
        Reset frames per integration
    frame_time : float
        Frame time (in seconds)
    frames_per_group : int, optional
        Frames per group -- always one except brown dwarves

    Returns
    -------
    duration_time : float
        Duration time (in seconds).
    """

    duration_time = frame_time * (num_groups * frames_per_group + num_reset_frames) * num_integrations

    return duration_time


def calc_exposure_time(num_integrations, ramp_time):
    """Calculates exposure time (or photon collection duration as told
    by APT.)

    Parameters
    ----------
    num_integrations : int
        Integrations per exposure.
    ramp_time : float
        Ramp time (in seconds).

    Returns
    -------
    exposure_time : float
        Exposure time (in seconds).
    """

    exposure_time = num_integrations * ramp_time

    return exposure_time


def calc_frame_time(num_columns, num_rows, num_amps, instrument):
    """Calculates the frame time for a given instrument/readmode/subarray.

    Parameters
    ----------
    num_columns : int
        Number of columns
    num_rows : int
        Number of rows
    num_amps : int
        Amplifiers reading data
    instrument : str
        The instrument

    Returns
    -------
    frame_time : float
        The frame time (in seconds)
    """

    num_columns, num_amps, num_rows = int(num_columns), int(num_amps), int(num_rows)

    if instrument == 'nirspec':
        n = 2
    if instrument in ['nircam', 'niriss']:
        n = 1

    frame_time = (num_columns / num_amps + 12) * (num_rows + n) * (1e-5)

    return frame_time


def calc_groups_from_exp_time(max_exptime_per_int, frame_time):
    """Given the maximum saturation time, calculates the number of
    frames per group.

    Parameters
    ----------
    max_exptime_per_int : float
        The maximum number of seconds an integration can last before
        it's oversaturated
    frame_time : float
        The time per frame

    Returns
    -------
    groups : int
        The required number of groups.
    """

    groups = np.floor(max_exptime_per_int / frame_time)

    return groups


def calc_integration_time(num_groups, frame_time, frames_per_group, num_skips):
    """Calculates the integration time.

    Parameters
    ----------
    num_groups : int
        Groups per integration.]
    frame_time : float
        Frame time (in seconds)
    frames_per_group : int
        Frames per group -- always 1 except maybe brown dwarves
    num_skips : int
        Skips per integration -- always 0 except maybe brown dwarves

    Returns
    -------
    integration_time : float
        Integration time (in seconds)
    """

    integration_time = (num_groups * (frames_per_group + num_skips) - num_skips) * frame_time

    return integration_time


def calc_num_integrations(transit_time, num_groups, num_reset_frames, frame_time, frames_per_group):
    """Calculates number of integrations required.

    Parameters
    ----------
    transit_time : float
        The time of the transit (in hours)
    num_groups : int
        Groups per integration
    num_reset_frames : int
        Number of reset frames per integration
    frame_time : float
        The frame time (in seconds)
    frames_per_group : int
        Frames per group -- always 1 except maybe for brown dwarves

    Returns
    -------
    num_integrations : float
        The required number of integraitions.
    """

    num_integrations = math.ceil((float(transit_time) * 3600) / (frame_time * (num_groups * frames_per_group + num_reset_frames)))

    return num_integrations


def calc_observation_efficiency(exposure_time, duration_time):
    """Calculates the observation efficiency.

    Parameters
    ----------
    exposure_time : float
        Exposure time (in seconds).
    duration_time : float
        Duration time (in seconds).

    Returns
    -------
    observation_efficiency : float
        Observation efficiency.
    """

    observation_efficiency = exposure_time / duration_time

    return observation_efficiency


def calc_ramp_time(integration_time, num_reset_frames, frame_time):
    """Calculates the ramp time -- or the integration time plus overhead for resets.

    Parameters
    ----------
    integration_time : float
        Integration time (in seconds)
    num_reset_frames : int
        Rest frames per integration
    frame_time : float
        Frame time (in seconds)

    Returns
    -------
    ramp_time : float
        Ramp time (in seconds).
    """

    ramp_time = integration_time + (num_reset_frames - 1) * frame_time

    return ramp_time


def convert_saturation(max_saturation, saturation_mode, instrument, infile, target_acq_mode=False):
    """Converts full well fraction to a saturation in counts OR
    provides the max fullwell for TA mode.

    Parameters
    ----------
    max_saturation : float
        Either a full well fraction or counts
    saturation_mode : str
        ``well`` or ``counts``
    instrument : str
        The instrument
    infile : str
        The path to the data file
    target_acq_mode : bool, optional
        Whether or not it's TA mode

    Returns
    -------
    max_saturation : float
        The fullwell to use in counts.
    """

    with open(infile) as f:
        data = json.load(f)

    instrument_dict = data['fullwell']

    if saturation_mode == 'well':
        max_saturation = float(max_saturation) * float(instrument_dict[instrument])

    if target_acq_mode:
        max_saturation = instrument_dict[instrument]

    return max_saturation


def interpolate_from_pandeia(magnitude, instrument, filt, subarray, model, band, frame_time, saturation_level, infile, target_acq_mode=False):
    """Interpolates the precalculated ``pandeia`` data to estimate the
    saturation limit.

    Parameters
    ----------
    magnitude : float
        The magnitude of the source. (Takes between 4.5-12.5)
    instrument : str
       The instrument, allowable ``miri``, ``niriss``, ``nirspec``,
       ``nircam``
    filt : str
        The filter
    subarray : str
        The subarray
    model : str
        Phoenix model key
    band : str
        Magnitude band
    frame_time : float
        Frame time
    saturation_level : float
        The maximum fullwell saturation we'll allow
    infile : str
        The data file to use
    target_acq_mode : bool, optional
        Whether or not we're running this for TA

    Returns
    -------
    num_groups : int
        The number of groups that won't oversaturate the detector
    max_sat : int
        The maximum saturation level reached by that number of groups
    """

    # Create the dictionaries for each filter and select out the prerun data
    with open(infile) as f:
        data = json.load(f)

    ta_or_sci = 'sci_sat'

    if target_acq_mode:
        ta_or_sci = 'ta_sat'

    # The data
    magnitudes = np.array(data['mags'])
    saturation = data[ta_or_sci][instrument][filt][subarray][model]
    log_saturation = np.log10(saturation)

    # Interpolate the given magnitude
    func_log = interpolate.interp1d(magnitudes, log_saturation)
    max_log_saturation = func_log(float(magnitude))
    max_saturation = 10**(max_log_saturation)

    # Figure out what it means in wake of the given saturation lvl
    max_exptime = float(saturation_level) / max_saturation

    # Calculate the nearest number of groups
    num_groups = calc_groups_from_exp_time(max_exptime, frame_time)

    # Can't have zero groups
    num_groups = num_groups or 1

    return num_groups, max_saturation


def map_to_ta_modes(instrument, max_num_groups, min_num_groups):
    """Turns the min/max groups into the closest allowable TA group
    mode.

    Parameters
    ----------
    instrument : str
        The instrument
    max_num_groups : int
        The maximum number of groups without oversaturating
    min_num_groups : int
        The groups needed to hit the target SNR

    Returns
    -------
    min_ta_groups : int
        The min possible groups to hit target SNR
    max_ta_groups : int
        The max possible groups before saturation
    """

    # Allowable group modes for each instrument
    groups = {'miri': [3, 5, 9, 15, 23, 33, 45, 59, 75, 93, 113, 135, 159, 185, 243, 275, 513],
              'niriss': [3, 5, 7, 9, 1, 13, 15, 17, 19],
              'nirspec': [3],
              'nircam': [3, 5, 9, 17, 33, 65]}

    # Match the literal min and max groups to the nearest mode.
    allowable_groups = groups[instrument]
    min_ta_groups = min(allowable_groups, key=lambda x: abs(x - min_num_groups))
    max_ta_groups = min(allowable_groups, key=lambda x: abs(x - max_num_groups))

    # Unless it was oversaturated from the get-go OR there aren't enough groups
    # for SNR
    if min_num_groups == 0:
        min_ta_groups = 0
        max_ta_groups = 0
    if min_num_groups > max(allowable_groups):
        min_ta_groups = -1
        max_ta_groups = 0

    return max_ta_groups, min_ta_groups


def min_num_groups_for_sat(magnitude, instrument, filt, subarray, model, band, infile):
    """Estimates the minimum number of groups to reach target acq
    saturation requirements.

    Parameters
    ----------
    magnitude : float
        Magnitude of star.
    instrument : str
        Instrument.
    filt : str
        Filter.
    subarray : str
        Subarray.
    model : str
        Phoenix model key.
    band : str, currently unused?
        The band -- right now only k sooo?
    infile : str
        The file with the pandeia data.

    Returns
    -------
    min_num_groups_for_sat : int
        The minimum number of groups to reach target snr.
    """

    with open(infile) as f:
        data = json.load(f)

    # Match to closest magnitude
    magnitudes = [float(i) for i in data['mags']]
    closest_magnitude = min(magnitudes, key=lambda x: abs(x - float(magnitude)))
    index = magnitudes.index(closest_magnitude)

    # Match to data
    minimum_num_groups_for_sat = data['ta_snr'][instrument][filt][subarray][model][index]

    return minimum_num_groups_for_sat


def perform_calculation(params, frames_per_group=1, num_skips=0):
    """Calculates all of the outputs and puts them in a dictionary for
    easy access.

    Parameters
    ----------
    params : dict
        Dictionary of all the needed parameters. Must include:
        ``obs_time``, ``n_group``, ``mag``, ``mod``, ``band``,
        ``filt``, ``filt_ta``, ``ins``, ``subarray``, ``subarray_ta``,
        ``sat_mode``, ``max_sat``, ``infile``
    frames_per_group : int, optional
        The number of frames -- almost always 1
    num_skips: int, optional
        Number of skips -- almost always 0

    Returns
    -------
    params : dict, str
        Dictionary of outputs and inputs. If the calculation throws an
        error it will return a string error message instead
    """

    # TARGET ACQ
    ta_frame_time = set_frame_time(
        params['infile'], params['ins'], params['subarray_ta'], target_acq_mode=True)

    max_saturation_ta_level = convert_saturation(
        params['sat_max'], params['sat_mode'], params['ins'], params['infile'], target_acq_mode=True)

    max_num_groups, saturation_rate_ta = interpolate_from_pandeia(
        params['mag'], params['ins'], params['filt_ta'], params['subarray_ta'], params['mod'], params['band'],
        ta_frame_time, max_saturation_ta_level, params['infile'], target_acq_mode=True)

    min_num_groups = min_num_groups_for_sat(
        params['mag'], params['ins'], params['filt_ta'], params['subarray_ta'], params['mod'], params['band'], params['infile'])

    min_ta_groups, max_ta_groups = map_to_ta_modes(params['ins'], max_num_groups, min_num_groups)

    duration_time_ta_min = calc_duration_time(min_ta_groups, 1, 1, ta_frame_time)

    duration_time_ta_max = calc_duration_time(max_ta_groups, 1, 1, ta_frame_time)

    # Science obs
    # Figure out the rows/cols/amps/pixel_size and saturation
    instrument_params = set_params_from_instrument(params['ins'], params['subarray'])
    frame_time = set_frame_time(params['infile'], params['ins'], params['subarray'])

    # Run all the calculations
    num_rows, num_columns, num_amps, pixel_size, frame_time, num_reset_frames = instrument_params
    params['sat_max'] = convert_saturation(params['sat_max'], params['sat_mode'], params['ins'], params['infile'])

    # Calculate countrate and n_groups if it isn't supplied
    num_groups, saturation_rate = interpolate_from_pandeia(
        params['mag'], params['ins'], params['filt'], params['subarray'], params['mod'], params['band'],
        frame_time, params['sat_max'], params['infile'])
    if str(params['n_group']) == 'optimize':
        params['n_group'] = int(num_groups)
    else:
        params['n_group'] = int(float(params['n_group']))

    # Aditional helpful params
    # Calculate times/ramps/etc
    integration_time = calc_integration_time(params['n_group'], frame_time, frames_per_group, num_skips)
    ramp_time = calc_ramp_time(integration_time, num_reset_frames, frame_time)

    # Calculate nubmer of integrations (THE MEAT)
    num_integrations = calc_num_integrations(params['obs_time'], params['n_group'], num_reset_frames, frame_time, frames_per_group)

    # Other things that may come in handy who knows?
    exposure_time = calc_exposure_time(num_integrations, ramp_time)
    duration_time = calc_duration_time(params['n_group'], num_integrations, num_reset_frames, frame_time, frames_per_group)
    observation_efficiency = calc_observation_efficiency(exposure_time, duration_time)

    # Update params with new information
    params['duration_time'] = round(duration_time / 3600, 3)
    params['duration_time_ta_max'] = duration_time_ta_max
    params['duration_time_ta_min'] = duration_time_ta_min
    params['exposure_time'] = round(exposure_time / 3600, 3)
    params['frames_per_group'] = frames_per_group
    params['frame_time'] = round(frame_time, 3)
    params['integration_time'] = round(integration_time, 3)
    params['max_saturation_prediction'] = round(saturation_rate * frame_time * params['n_group'], 3)
    params['max_saturation_ta'] = round(saturation_rate_ta * ta_frame_time * max_ta_groups, 3)
    params['min_saturation_ta'] = round(saturation_rate_ta * ta_frame_time * min_ta_groups, 3)
    params['max_ta_groups'] = int(max_ta_groups)
    params['min_ta_groups'] = int(min_ta_groups)
    params['num_amps'] = num_amps
    params['num_columns'] = num_columns
    params['num_integrations'] = num_integrations
    params['num_reset_frames'] = num_reset_frames
    params['num_rows'] = num_rows
    params['num_skips'] = num_skips
    params['observation_efficiency'] = round(observation_efficiency, 3)
    params['ramp_time'] = round(ramp_time, 3)
    params['ta_frame_time'] = ta_frame_time

    return params


def set_params_from_instrument(instrument, subarray):
    """Sets/collects the running parameters from the instrument.

    Parameters
    ----------
    instrument : str
        Instrument, options are ``nircam``, ``niriss``, ``nirpec``, and
        ``miri``
    subarray : str
        Subarray mode

    Returns
    -------
    rows : int
        The number of pixels per row.
    cols : int
        The number of columns per row.
    amps : int
        The number of amplifiers.
    pixel_size : int
        The pixel size.
    frame_time : float
        The frame time.
    num_reset_frames : int
        The number of reset frames.
    """

    num_reset_frames = 1

    if instrument == 'nirspec':

        pixel_size = (40e-4)**2
        if subarray == 'sub2048':
            rows, cols = 2048, 32
        elif subarray in ['sub1024a', 'sub1024b']:
            rows, cols = 1024, 32
        elif subarray == 'sub512':
            rows, cols = 512, 32
        amps = 1  # 4 if not NRSRAPID????
        frame_time = calc_frame_time(cols, rows, amps, instrument)

    elif instrument == 'nircam':

        pixel_size = (18e-4)**2
        if subarray == 'full':
            rows, cols, amps = 2048, 2048, 4
        elif subarray == 'subgrism256':
            rows, cols, amps = 256, 256, 1
        elif subarray == 'subgrism128':
            rows, cols, amps = 128, 2048, 1
        elif subarray == 'subgrism64':
            rows, cols, amps = 64, 2048, 1
        frame_time = calc_frame_time(cols, rows, amps, instrument)

    elif instrument == 'miri':

        pixel_size = (25e-4)**2
        if subarray == 'slitlessprism':
            rows, cols, frame_time = 416, 72, .159
        amps = 4
        num_reset_frames = 0

    elif instrument == 'niriss':
        pixel_size = (40e-4)**2
        if subarray == 'substrip96':
            rows, cols = 2048, 96
        elif subarray == 'substrip256':
            rows, cols = 2048, 256
        amps = 1
        frame_time = calc_frame_time(cols, rows, amps, instrument)

    return rows, cols, amps, pixel_size, frame_time, num_reset_frames


def set_frame_time(infile, instrument, subarray, target_acq_mode=False):
    """Assign the appropriate frame time based on the instrument and
    subarray. For now, modes are implied.

    Parameters
    ----------
    infile: str
        The path to the data file.
    instrument : str
        The instrument : ``miri``, ``niriss``, ``nirspec``, or
        ``nircam``
    subarray : str
        The subarray
    target_acq_mode : bool
        Whether this is for TA or not.

    Returns
    -------
    frame_time : float
        The frame time for this instrument/subarray combo.
    """

    # Read in dict with frame times
    with open(infile) as f:
        frame_time = json.load(f)['frame_time']

    if target_acq_mode:
        frame_time = frame_time[instrument]['ta'][subarray]
    else:
        frame_time = frame_time[instrument][subarray]

    return frame_time
