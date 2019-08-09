import numpy as np

from exoctk.utils import get_target_data

def deg2HMS(ra='', dec='', round=True):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = np.round(((abs((dec-deg)*60)-decM)*60), 2)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1}:{2}:{3}'.format(ds, deg, decM, decS)

  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = np.round(((((ra/15)-raH)*60)-raM)*60, 2)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1}:{2}:{3}'.format(rs, raH, raM, raS)

  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

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

    ra_hms, dec_hms = deg2HMS(ra=ra_deg, dec=dec_deg)

    return ra_hms, dec_hms
