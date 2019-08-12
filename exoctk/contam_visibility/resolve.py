import numpy as np
from astropy.coordinates import SkyCoord
from exoctk.utils import get_target_data

def resolve_target(targetName):
    """
    Gets the RA and Dec of your target.

    Parameters
    ----------
    targetName : str
        The name of your target.

    Returns
    -------
    ra : str
        The Right Ascension (RA) of your target in HH:MM:SS.
    dec : str
        The Declination (Dec) of your target in HH:MM:SS.
    """
    data, url = get_target_data(targetName)
    ra_deg = data['RA']
    dec_deg = data['DEC']

    sc = SkyCoord(ra_deg, dec_deg, unit='deg')
    ra_hms = '{}:{}:{}'.format(int(sc.ra.hms[0]), \
                               int(sc.ra.hms[1]), \
                               np.round(sc.ra.hms[2], 2))
    dec_hms = '{}:{}:{}'.format(int(sc.dec.dms[0]), \
                                int(sc.dec.dms[1]), \
                                np.round(sc.dec.dms[2], 2))

    return ra_hms, dec_hms
