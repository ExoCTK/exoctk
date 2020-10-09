import numpy as np

from astropy.coordinates import SkyCoord
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show

from exoctk.contam_visibility import visibilityPA as vpa
from exoctk.contam_visibility import field_simulator as fs
from exoctk.contam_visibility import contamination_figure as cf

def main():
    """ Wrapper to the field simulator and contamination figure generator.
    """

    # User inputs
    ra = input('Please input the Right Ascension of your target in decimal degrees (The more decimal places the better) : \n')
    # Making sure RA input is correct before we continue
    if (float(ra) < 0) or (float(ra) > 360):
        print('RA should be between 0 and 360 decimal degrees. Got {}. Starting over...'.format(ra))
        main()

    dec = input('Please input the Declination of your target in decimal degrees (The more decimal places the better) : \n')
    # Making sure DEC input is correct
    if (float(dec) < -90) or (float(dec) > 90):
        print('DEC should be between -90 and +90 in decimal degrees. Got {}. Starting over...'.format(dec))
        main()

    companion = input('Any companion not in IRSA`s 2MASS Point-Source Catalog that should be considered for contamination? If no, just press Enter. If yes, please enter the following (comma-separated, no spaces): \n RA offset ("), DEC offset ("), 2MASS J (mag), H (mag) and Ks (mag) \n ').split(',')
    print(len(companion))
    if len(companion)==5:
        binComp = [int(param) for param in companion]
    elif len(companion)==1:
        binComp = ''
    elif (len(companion) < 5) & (len(companion) > 1):
        print('Companion information is incomplete. Starting over...')
        main()

    instrument = input('Please input the instrument your target will be observed with (Case-sensitive. Can be: NIRISS, MIRI, NIRCam F322W2, or NIRCam F444W) : \n')
    # Making sure the instrument input is correct before we continue
    possible_instruments = ['NIRISS', 'MIRI', 'NIRCam F322W2', 'NIRCam F444W']
    if instrument not in possible_instruments:
        print('It looks like your last input (instrument) had a typo. The instrument input is case-sensitive and can only be NIRISS, MIRI, NIRCam F322W2, or NIRCam F444W. Starting over...')
        main()

    # Getting the bad PAs from visibility calculator for shading purposes
    instrument_vpa = instrument.split(' ')[0]
    paMin, paMax, gd, fig, table, grouped_badPAs = vpa.using_gtvt(ra, dec, instrument_vpa)

    # Converting RA, DEC from decimal degrees (ExoMAST) to HMSDMS
    sc = SkyCoord(ra, dec, unit='deg')
    ra_dec = sc.to_string('hmsdms')
    ra_hms, dec_dms = ra_dec.split(' ')[0], ra_dec.split(' ')[1]

    # Generating cube with a field for every Aperture Position Angle (APA)
    cube = fs.fieldSim(ra_hms, dec_dms, instrument, binComp)

    title_ra, title_dec = str(np.round(float(ra), 3)), str(np.round(float(dec), 3))
    # Generating Bokeh figure `fig` that plots contamination levels at every APA
    plot = cf.contam(cube, instrument, targetName=' {}, {} (RA, DEC)'.format(title_ra, title_dec), badPAs=grouped_badPAs)

    show(plot)

if __name__ == "__main__":
    main()
