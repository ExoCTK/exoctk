from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad


def resolve_target(target_name):
    """Resolve a target to its SIMBAD ICRS/J2000 coordinates.

    Parameters
    ----------
    target_name : str
        Target identifier accepted by SIMBAD's basic identifier query.

    Returns
    -------
    tuple of float
        Right ascension and declination in decimal degrees.

    Raises
    ------
    ValueError
        If SIMBAD cannot resolve ``target_name``.
    """

    result = Simbad.query_object(target_name)
    if result is None or len(result) == 0:
        raise ValueError(f"SIMBAD could not resolve '{target_name}'.")
    coordinate = SkyCoord.guess_from_table(result)
    if not coordinate.isscalar:
        coordinate = coordinate[0]

    return float(coordinate.ra.deg), float(coordinate.dec.deg)
