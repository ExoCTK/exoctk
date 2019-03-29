from exoctk.utils import get_target_data

def resolve_target(targetName):
    data = get_target_data(targetName)
    ra = data['RA']
    dec = data['DEC']

    return ra, dec
