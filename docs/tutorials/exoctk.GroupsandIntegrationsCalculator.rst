.. _GroupsandIntegrationsCalculator:

Groups & Integrations Calculator
================================

.. contents::

This is a quick demo for the Groups and Integrations Calculator (i.e. the `exoctk.groups_integrations` subpackage). This demo can run with the `exoctk-3.8` `conda` environment and a proper installation of the `exoctk` software package and the `EXOCTK_DATA` data package.  For instructions on how to download and install these things, see the `installation instructions <https://github.com/ExoCTK/exoctk#installation>`_.

The Groups and Integrations Calculator was written as an exoplanet-focused alternative to the JWST ETC (also known as `pandeia`). It's power comes from:

1. Using pre-sampled data and interpolating for a speedy calculation.
2. Having the power to provide you with the optimal configuration for your observation -- instead of just letting you guess and check configurations.

This notebook has a few sections for clarity:

- Imports and Setup
- Input Parameter Space
- Building an Input Dictionary
- Running the Calculator
- Exploring the Outputs
- Two Examples for Batch Runs


Imports and Setup
-----------------

In addition to installing `exoctk`, you'll need the pre-sampled data that the groups and integrations calculator runs on, in a `json` format.  If you followed the installation instructions linked at the top of this notebook, you should have this under the environment variable `EXOCTK_DATA`, and we'll just take it from your `path` -- otherwise you'll need to specify the path in this cell.

Also -- **DON'T WORRY** if you get a big scary warning about the `matplotlib` backend. We need that in there to run server side when we host the website -- but it won't affect any of this demo of your general `exoctk` usage. If it's truly too hideous run this cell twice and warning should vanish.

.. code:: ipython3

    # Check if you have the environment variable
    import os
    EXOCTK_DATA = os.environ.get('EXOCTK_DATA')
    if EXOCTK_DATA == '':
        GROUPS_INTEGRATIONS_DIR = '/path/to/your/local/copy/integrations_groups_data.json'
    else:
        GROUPS_INTEGRATIONS_DIR = GROUPS_INTEGRATIONS_DIR = os.path.join(EXOCTK_DATA, 'groups_integrations/groups_integrations.json')

.. code:: ipython3

    # Imports
    import json
    from exoctk.groups_integrations.groups_integrations import perform_calculation


Input Parameter Space
---------------------

We are often adding functionality and modes to this -- so instead of providing a static list of what models, instrument modes, etc., are presently feasible, here's a little demo on how to walk through the data file yourself and see what parameters are possible.

.. code:: ipython3

    # Open the file
    with open(GROUPS_INTEGRATIONS_DIR) as f:
        parameter_space = json.load(f)

    # Print the TA and Science observation modes
    modes = [('science', 'sci_sat'), ('TA', 'ta_sat')]
    for mode in modes:
        print('\n')
        print('Available modes for {}:\n'.format(mode[0]))
        instruments = list(parameter_space[mode[1]].keys())
        print('Instruments: {}\n'.format(instruments))
        for instrument in instruments:
            filters = list(parameter_space[mode[1]][instrument].keys())
            print('Filters for {}: {}'.format(instrument, filters))
            subarrays = list(parameter_space[mode[1]][instrument][filters[0]].keys())
            print('Subarrays for {}: {}\n'.format(instrument, subarrays))

    # Print the available Phoenix models
    print('\n')
    print('Phoenix model keys :')
    print('--------------------')
    print(list(parameter_space['ta_sat'][instruments[-1]][filters[0]][subarrays[0]].keys()))

    print('\n')
    print('Magnitude sampling :')
    print('--------------------')
    print(parameter_space['mags'])

.. code:: txt

    Available modes for science:

    Instruments: ['nirspec', 'niriss', 'miri', 'nircam']

    Filters for nirspec: ['f070lp_g140m', 'f070lp_g140h', 'f170lp_g235h', 'f170lp_g235m', 'f290lp_g395m', 'f100lp_g140h', 'f290lp_g395h', 'f100lp_g140m', 'clear_prism']
    Subarrays for nirspec: ['sub512', 'sub2048', 'sub1024a', 'sub1024b']

    Filters for niriss: ['soss']
    Subarrays for niriss: ['substrip256', 'substrip96']

    Filters for miri: ['lrs']
    Subarrays for miri: ['slitlessprism']

    Filters for nircam: ['f322w2', 'f444w', 'f277w', 'f356w']
    Subarrays for nircam: ['subgrism128', 'subgrism64', 'full', 'subgrism256']


    Available modes for TA:

    Instruments: ['nirspec', 'niriss', 'miri', 'nircam']

    Filters for nirspec: ['f140x', 'f110w', 'clear']
    Subarrays for nirspec: ['full', 'sub32', 'sub2048']

    Filters for niriss: ['f480m']
    Subarrays for niriss: ['subtasoss']

    Filters for miri: ['f560w', 'f100w', 'f1500w']
    Subarrays for miri: ['slitlessprism']

    Filters for nircam: ['f335m']
    Subarrays for nircam: ['sub32tats']


    Phoenix model keys :
    --------------------
    ['g5v', 'a1v', 'a3v', 'k5iii', 'm0i', 'm0v', 'g5iii', 'f5v', 'k0v', 'g5i', 'k2v', 'frame_time', 'm2i', 'k0i', 'g0iii', 'f5i', 'm0iii', 'a0v', 'g8v', 'a0i', 'mag', 'k7v', 'g2v', 'm5v', 'g0i', 'a5i', 'k5v', 'f8v', 'f0i', 'f0v', 'f2v', 'k5i', 'a5v', 'k0iii', 'g0v', 'm2v']


    Magnitude sampling :
    --------------------
    [4.5, 6.5, 8.5, 10.5, 12.5]


Building an Input Dictionary
----------------------------

Running the groups and integrations calculator requires a dictionary of inputs. This section will go through an example input dictionary and what the limits on the parameters are.

.. code:: ipython3

    # Initialize the dictionary
    parameters = {}

    # Source parameters
    parameters['mag'] = 10 # 4.5 <= float <= 12.5
    parameters['band'] = 'k' # only K band vega mag for now
    parameters['mod'] = 'g2v' # Phoenix model per last section

    # Observation specifics
    parameters['obs_time'] = 5 # positive float, in hours
    parameters['n_group'] = 'optimize' # 'optimize', or positive integer

    # Detector setup -- within the modes of the last section
    parameters['ins'] = 'nircam'
    # For science observation
    parameters['filt'] = 'f444w'
    parameters['subarray'] = 'subgrism256'
    # And target acquisition
    parameters['filt_ta'] = 'f335m'
    parameters['subarray_ta'] = 'sub32tats'

    # Saturation level
    parameters['sat_mode'] = 'well' # 'well', for full well fraction, or 'counts'
    parameters['sat_max'] = .95 # < 1 for fullwell, and a positive integer always

    # And feed in the data file
    parameters['infile'] = GROUPS_INTEGRATIONS_DIR
    input_dict = parameters.copy()


Running the Calculator
----------------------

Now, running the calculator is relatively straightforward. We leaned on the `pandeia` function convention -- so feeding our inputs into `perform_calculation` returns a dictionary of input parameters (a stripped down `pandiea` scene) as well as the results of the calculation. (Note that `perform_calculation` updates the `parameters` dictionary -- so that object will be changed once you run the function.)

.. code:: ipython3

    # Bookeeping for new/old parameters
    inputs = list(parameters.keys())

    # Perform the calculation
    results = perform_calculation(parameters)
    for key in results:
        if key in inputs:
            if key == 'infile':
                # hackers
                print('The input of infile was REDACTED!')
            else:
                print('The input of {} was {}.'.format(key, results[key]))
        else:
            print('The result of {} was {}'.format(key, results[key]))

.. code:: txt

    The input of mag was 10.
    The input of band was k.
    The input of mod was g2v.
    The input of obs_time was 5.
    The input of n_group was 149.
    The input of ins was nircam.
    The input of filt was f444w.
    The input of subarray was subgrism256.
    The input of filt_ta was f335m.
    The input of subarray_ta was sub32tats.
    The input of sat_mode was well.
    The input of sat_max was 55195.0.
    The input of infile was REDACTED!
    The result of n_col was 256
    The result of n_row was 256
    The result of n_amp was 1
    The result of n_reset was 1
    The result of n_frame was 1
    The result of n_skip was 0
    The result of t_frame was 1.347
    The result of t_int was 200.657
    The result of t_ramp was 200.657
    The result of n_int was 90
    The result of t_exp was 5.016
    The result of t_duration was 5.05
    The result of obs_eff was 0.993
    The result of ta_t_frame was 0.01496
    The result of min_ta_groups was 33
    The result of max_ta_groups was 3
    The result of t_duration_ta_min was 0.50864
    The result of t_duration_ta_max was 0.05984
    The result of max_sat_prediction was 55171.429
    The result of max_sat_ta was 4613.154
    The result of min_sat_ta was 50744.699


Exploring the Outputs
---------------------

If you aren't quite familiar with the intricacies of a JWST observation, we'll unpack these results briefly.

Every JWST observation has a number of groups and integrations. Groups are how many frames you can pack into an integration, and generally this is goverened by how quickly your target will saturate on the detector. With every integration, there is overhead time added into the observation, and less time observing your actual transit. So, once the calculator has figured out the maximum possible groups before the detector is saturated, it will return that number of groups, how many integrations of that pattern it takes to fill up your whole transit time, and some additional helpful parameters like what's the maximum saturation you might reach during this observation, how long the actual scheme will take (since there will be a slightly different rounding when everything is added up), how efficient your observation is, etc.

For target acquisition, the efficiency is less in question than the ability of the detector to hit the minimum SNR required. This provides a reccomendation, so you know the minimum and maxmimum possible groups you can use, and can make an informed decision.

.. code:: ipython3

    # So let's make a nice printed summary
    print('The total time for science + TA observation scheme is {}-{} hours.'.format(
        results['t_duration']+results['t_duration_ta_max'], results['t_duration']+results['t_duration_ta_min']))
    print('You need {} groups and {} integrations for the science observation.'.format(
        results['n_group'], results['n_int']))
    print('You need between {} and {} groups for target acquisition.'.format(
        results['max_ta_groups'], results['min_ta_groups']))
    print('We estimate your science observation will reach at most {} counts -- how close were we to your cutoff of {}?'.format(
        results['max_sat_prediction'], results['sat_max']))
    print('With this observation scheme {}% of the observation will be science data.'.format(results['obs_eff']*100))

.. code:: txt

    The total time for science + TA observation scheme is 5.10984-5.55864 hours.
    You need 149 groups and 90 integrations for the science observation.
    You need between 3 and 33 groups for target acquisition.
    We estimate your science observation will reach at most 55171.429 counts -- how close were we to your cutoff of 55195.0?
    With this observation scheme 99.3% of the observation will be science data.


Two Examples for Batch Runs
---------------------------

So far we've shown you how to run this one off -- just like you would in the `web tool <https://exoctk.stsci.edu/groups_integrations>`_.  Here are two examples for running many calculations. Because the calculator is so light and only has a single parameter, it won't be particularly computationally expensive or logistically difficult to parallelize.

- Checking many potential observations in *one* instrument/mode/filter.
- Checking one transit in *many* instruments/modes/filters.


Checking many potential observations in *one* instrument/mode/filter
--------------------------------------------------------------------

.. code:: ipython3

    # Some imports/set up for ease of this part of the demo
    from multiprocessing import Pool
    p = Pool(4) # Feel free to set it higher but this will place nice with most laptops.

    from astropy.io import ascii
    import numpy as np

.. code:: ipython3

    # Say you have a table of sources you want to read in.
    # sources = ascii.read('path/to/source_table.csv')
    # Since we don't we'll just make one up with reasonable transit objects
    sources = {'mod': ['k0v', 'k5v', 'g5v', 'f5i', 'g0iii', 'f0i', 'k0iii', 'g2v', 'm0v', 'k5iii', 'm0i', 'g0v', 'g8v', 'f0v', 'g0i'],
               'obs_time': [3 + n*.15 for n in range(15)],
               'mag': [9 + n*.15 for n in range(15)]}

    # Now use this to create input dictionaries
    input_sources = []
    for index, elem in enumerate(sources['mod']):
        input_source = input_dict.copy()
        input_source['mod'] = elem
        input_source['obs_time'] = sources['obs_time'][index]
        input_source['mag'] = sources['mag'][index]
        input_sources.append(input_source)

    # And run it in parallel
    single_mode_results = p.map(perform_calculation, input_sources)

.. code:: ipython3

    # And explore the output
    obs_eff = [result['obs_eff'] for result in single_mode_results]
    indeces = np.where(obs_eff == np.max(obs_eff))[0]
    bests = [single_mode_results[index] for index in indeces]
    for best in bests:
        print('One of the best sources is {}, {}, {} at {} efficiency.'.format(best['mod'], best['obs_time'], best['mag'], best['obs_eff']))
        print('(This means {} groups and {} integrations.)'.format(best['n_group'], best['n_int']))

.. code:: txt

    One of the best sources is g0i, 5.1, 11.1 at 0.998 efficiency.
    (This means 406 groups and 34 integrations.)


Checking one transit in *many* instruments/modes/filters
--------------------------------------------------------

.. code:: ipython3

    # Let's take the LEAST efficient observation from the last example
    # We'll see if it plays more nicely with another instrument, filter, or subarray.
    worst = single_mode_results[np.where(obs_eff == np.min(obs_eff))[0][0]]
    print("We're starting at a baseline of {} efficiency.".format(worst['obs_eff']))
    print('Source : {}, {}, {} mode.'.format(worst['mod'], worst['mag'], worst['obs_time']))
    print('Mode : {}, {}, {}'.format(worst['ins'], worst['filt'], worst['subarray']))
    for key in ['mod', 'mag', 'obs_time']:
        input_dict[key] = worst[key]

    # We'll call back to the parameter_space dictionary we walked through to look at available modes
    modes = []
    for ins in parameter_space['sci_sat'].keys():
        for filt in parameter_space['sci_sat'][ins].keys():
            for sub in parameter_space['sci_sat'][ins][filt].keys():
                input_mode = input_dict.copy()

                # Alter the science setup
                input_mode['ins'] = ins
                input_mode['filt'] = filt
                input_mode['subarray'] = sub

                # And we care less about TA so pick the default for each instrument
                input_mode['filt_ta'] = list(parameter_space['ta_sat'][ins].keys())[0]
                input_mode['subarray_ta'] = list(parameter_space['ta_sat'][ins][input_mode['filt_ta']].keys())[0]
                modes.append(input_mode)

    # And run it in parallel
    single_source_results = p.map(perform_calculation, modes)

.. code:: txt

    We're starting at a baseline of 0.984 efficiency.
    Source : k0v, 9.0, 3.0 mode.
    Mode : nircam, f444w, subgrism256

.. code:: ipython3

    # And, again explore the output
    obs_eff = [result['obs_eff'] for result in single_source_results]
    indeces = np.where(obs_eff == np.max(obs_eff))[0]
    bests = [single_source_results[index] for index in indeces]
    for best in bests:
        print('The best mode is {}, {}, {} at {} efficiency.'.format(best['ins'], best['filt'], best['subarray'], best['obs_eff']))
        print('(This means {} groups and {} integrations.)'.format(best['n_group'], best['n_int']))

.. code:: txt

    The best mode is nircam, f444w, subgrism64 at 0.996 efficiency.
    (This means 234 groups and 135 integrations.)