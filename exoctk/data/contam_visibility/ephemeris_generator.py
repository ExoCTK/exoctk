"""
Script to generate ephemeris for contamination & visibility tool

>>> from exoctk.data.contam_visibility.ephemeris_generator import generate_ephemeris
>>> ephemeris_generator('2021-12-26', '2022-10-31', 'JWST_ephemeris_2021-2022.txt')
"""
from datetime import datetime
import numpy as np
import os
import pandas as pd
import requests

DEFAULT_OUTNAME_PATH = os.path.dirname(os.path.abspath(__file__))
URL = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND=-170&OBJ_DATA=%27NO%27&EPHEM_TYPE=VECTORS&START_TIME=%27{}%27&STOP_TIME=%27{}%27&CENTER=%27500@10%27&STEP_SIZE=%271%20DAYS%27&CSV_FORMAT=%27YES%27&VEC_TABLE=%272%27&REF_SYSTEM=%27ICRF%27&REF_PLANE=%27FRAME%27&VEC_CORR=%27LT%27&OUT_UNITS=%27KM-S%27&VEC_LABELS=%27YES%27&VEC_DELTA_T=%27NO%27'


def convert_ephemeris_to_df(ephemeris):

    start_index = np.where(ephemeris == '$$SOE')[0][0] + 1
    end_index = np.where(ephemeris == '$$EOE')[0][0]

    row_data = [row_data.split(',') for row_data in ephemeris[start_index:end_index]]
    result = [filter(None, row) for row in row_data]

    df = pd.DataFrame(result, columns=['JDTDB','Calendar Date (TDB)','X','Y','Z','VX', 'VY', 'VZ'])
    convert_dict = {'JDTDB': float,
                    'X': float, 'Y': float, 'Z': float,
                    'VX': float, 'VY': float, 'VZ': float}

    df = df.astype(convert_dict)

    return df


def get_ephemeris_data(start_date, end_date):
    """dates must be in format YYYY-MM-DD, returns text of ephemeris
    """ 
    try:
        url = URL.format(start_date, end_date)  # Get Horizons url for JWST ephemeris and add user specified dates
        eph_request = requests.get(url)
        ephemeris = np.array(eph_request.text.splitlines())
    except requests.exceptions.ConnectionError as err:
        print('Issues connecting to JPL Horizons, check network connection.')

    return ephemeris


def generate_ephemeris(start_date, end_date, outname, path=None):
    """Read and write ephmeris in format for exoctk contam tool

    Parameters
    ----------
    eph : str
        path to ephemeris. Expecting files with names=['Time', 'x (km)', 'y (km)', 'z (km)', 'Dx (km/s)', 'Dy (km/s)', 'Dz (km/s)'].
        These columns map to [Calendar Date (TDB), X, Y, Z, VX, VY, VZ] JPL Horizons ephemeris columns. The ephemeris we are using to
        generate the "new" ephemeris is the data copied from JPL Horizons and the columns names changes to the expected names.
    outname : str
        Name to write ephemeris out to.
    path : str
        Path to write file out to, if one is not provided default is exoctk/data/contam_visibility/
    """

    eph = get_ephemeris_data(start_date, end_date)
    horizons_df = convert_ephemeris_to_df(eph)
    
    df = pd.DataFrame()
    df['Time'] = horizons_df['Calendar Date (TDB)']
    df['x (km)'] = horizons_df['X']
    df['y (km)'] = horizons_df['Y']
    df['z (km)'] = horizons_df['Z']
    df['Dx (km/s)'] = horizons_df['VX']
    df['Dy (km/s)'] = horizons_df['VY']
    df['Dz (km/s)'] = horizons_df['VZ']

    # Calculate vector magnitudes for distance and velocity
    df['Distance (km)'] = np.sqrt(df['x (km)'] **2 +  df['y (km)']**2 + df['z (km)']**2)
    df['Vel. (km/s)'] = np.sqrt(df['Dx (km/s)']**2 + df['Dy (km/s)']**2 + df['Dz (km/s)']**2)

    # Reformat Time
    df['Time'] = [reformat_time(date_str) for date_str in df['Time']]

    if path:
        df.to_csv(os.path.join(path, outname), index=False)
    else:
        df.to_csv(os.path.join(DEFAULT_OUTNAME_PATH, outname), index=False)


def reformat_time(date_str):
    """
    Parameters
    ----------
    data_str : str
        date string in format of "%Y-%b-%d %H:%M:%S.%f"
    """

    date_str = date_str.replace('A.D.', '')
    date_str = date_str.replace(' ', '')
    time_format = "%Y-%b-%d%H:%M:%S.%f"
    day = datetime.strptime(date_str, time_format)
    year = day.year
    day_of_year = day.strftime('%j')
    hour_minutes_seconds = day.strftime('%M:%S.%f')
    new_time = "{}-{}/{}".format(year, day_of_year, hour_minutes_seconds)

    return new_time