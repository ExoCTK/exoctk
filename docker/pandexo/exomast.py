import itertools
import os
import re
import urllib

import tornado.escape
from tornado.httpclient import AsyncHTTPClient


def build_target_url(target_name):
    '''Build restful api url based on target name.
    Parameters
        ----------
        target_name : string
            The name of the target transit.
        Returns
        -------
        target_url : string
    '''
    # Encode the target name string.
    encode_target_name = urllib.parse.quote(target_name, encoding='utf-8')
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/{}/properties/".format(encode_target_name)

    return target_url

async def get_canonical_name(target_name):
    '''Get ExoMAST prefered name for exoplanet.
        Parameters
        ----------
        target_name : string
            The name of the target transit.
        Returns
        -------
        canonical_name : string
    '''

    target_url = f"https://exo.mast.stsci.edu/api/v0.1/exoplanets/identifiers/?name={target_name}"

    http_client = AsyncHTTPClient()
    response = await http_client.fetch(url)
    planetnames = tornado.escape.json_decode(response.body)
    canonical_name = planetnames['canonicalName']

    return canonical_name

async def get_target_data(target_name):
    """
    Send request to exomast restful api for target information.
    Parameters
    ----------
    target_name : string
        The name of the target transit
    Returns
    -------
    target_data: json:
        json object with target data.
    """

    canonical_name = get_canonical_name(target_name)

    target_url = build_target_url(canonical_name)

    http_client = AsyncHTTPClient()
    response = await http_client.fetch(url)

    if r.status_code == 200:
        target_data = tornado.escape.json_decode(response.body)
    else:
        raise Exception('Whoops, no data for this target!')

    # Some targets have multiple catalogs
    # nexsci is the first choice.
    if len(target_data) > 1:
        # Get catalog names from exomast and make then the keys of a dictionary
        # and the values are its position in the json object.
        catalog_dict = {data['catalog_name']: index for index, data in enumerate(target_data)}

        # Parse based on catalog accuracy.
        if 'nexsci' in list(catalog_dict.keys()):
            target_data = target_data[catalog_dict['nexsci']]
        elif 'exoplanets.org' in list(catalog_dict.keys()):
            target_data = target_data[catalog_dict['exoplanets.org']]
        else:
            target_data = target_data[0]
    else:
        target_data = target_data[0]

    # Strip spaces and non numeric or alphabetic characters and combine.
    url = 'https://exo.mast.stsci.edu/exomast_planet.html?planet={}'.format(re.sub(r'\W+', '', canonical_name))

    return target_data, url
