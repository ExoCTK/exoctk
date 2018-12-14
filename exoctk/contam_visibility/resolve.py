from exoctk.utils import get_target_data

def resolve_target(targetName):
    data = get_target_data(targetName)
    ra = data['ra']
    dec = data['dec']

    return ra, dec

# from astroquery.simbad import Simbad


# def resolve_target(targetName):
#     """Resolve a target in SIMBAD

#     Parameters
#     ----------
#     targetName: str
#         The name of the target.

#     Returns
#     -------
#     tuple
#         The ra and dec of the target.
#     """
#     try:
#         target_info = Simbad.query_object(targetName)
#         targetRA = target_info['RA']
#         targetDEC = target_info['DEC']
#         ra = (targetRA[0].replace(' ', ':'))
#         dec = (targetDEC[0].replace(' ', ':'))
#         return ra, dec

#     except IndexError:
#         ra = 'unresolved'
#         dec = 'unresolved'
#         return ra, dec
