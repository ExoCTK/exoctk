from astroquery.simbad import Simbad


def resolve_target(targetName):
    """Resolve a target in SIMBAD

    Parameters
    ----------
    targetName: str
        The name of the target.

    Returns
    -------
    tuple
        The ra and dec of the target.
    """
    try:
        target_info = Simbad.query_object(targetName)
        targetRA = target_info['RA']
        targetDEC = target_info['DEC']
        ra = (targetRA[0].replace(' ', ':'))
        dec = (targetDEC[0].replace(' ', ':'))
        return ra, dec

    except IndexError:
        ra = 'unresolved'
        dec = 'unresolved'
        return ra, dec
