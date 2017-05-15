
from . import _exotransmit_wrapper

import inspect
import errno
import os
import requests
import sys

if sys.version_info.major >= 3:
    from urllib.parse import urljoin
else:
    from urlparse import urljoin

def exotransmit(**kwargs):
    """
    Run exotransmit.  The function can take any arguments for 
    `create_user_input`, `create_chem_selection` and `create_other_input` for 
    parameter names and descriptions.

    """

    check_user_input(**kwargs)
    create_user_input(**kwargs)
    create_chem_selection(**kwargs)
    check_other_input(**kwargs)
    create_other_input(**kwargs)

    _exotransmit_wrapper.exotransmit()

def check_user_input(base_dir=None, T_P_file='/T_P/t_p_800K.dat',
                      EOS_file='/EOS/eos_1Xsolar_cond.dat',
                      output_file='/Spectra/default.dat', g=9.8,
                      R_planet=6.4e+6, R_star=7.0e+8, P_cloud=0.0,
                      Rayleigh=1.0, **kwargs):

    if not base_dir:
        base_dir = os.path.abspath(os.curdir)

    # Need to catch missing files or else C will exit an interactive session
    if not os.path.isfile(base_dir + EOS_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                base_dir + EOS_file)

    if not os.path.isfile(base_dir + T_P_file):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                base_dir + T_P_file)

    if not os.path.isdir(base_dir + os.path.dirname(output_file)):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                base_dir + os.path.dirname(output_file))

    # other parameters need to be able to be cast to float
    float(g)
    float(R_planet)
    float(R_star)
    float(P_cloud)
    float(Rayleigh)

def create_user_input(base_dir=None, T_P_file='/T_P/t_p_800K.dat',
                      EOS_file='/EOS/eos_1Xsolar_cond.dat',
                      output_file='/Spectra/default.dat', g=9.8,
                      R_planet=6.4e+6, R_star=7.0e+8, P_cloud=0.0,
                      Rayleigh=1.0, **kwargs):
    """
    Create the userInput.in file needed to run Exo_Transmit.
    
    Parameters
    ----------
    base_dir: str
        The base working directory 
    T_P_file: str
        The location of the T-P profile file relative to the `base_dir`
    EOS_file: str
        The location of the Equation of State file relative to the `base_dir`
    output_file: str
        The location to save the output spectrum relaitve to the `base_dir`
    g: float
        The planet surface gravity in m / s^2
    R_planet: float
        The radius of the planet in m at the base of the atmosphere
    R_star: float
        The radius of the star in m 
    P_cloud: float
        The pressure at the top of the cloud in Pa, leave at zero for no
        cloud calculations
    Rayleigh: float
        Rayleigh scattering augmentation factor

    """

    if not base_dir:
        base_dir = os.path.abspath(os.curdir)

    with open(os.path.join(os.path.dirname(__file__),
                           'include/userInput.in')) as f:
        userinput_template = f.read()

    contents = userinput_template.format(base_dir, T_P_file, EOS_file,
                                         output_file, g, R_planet, R_star,
                                         P_cloud, Rayleigh)
    with open('userInput.in', 'w') as f:
        f.write(contents)

def create_chem_selection(CH4=True, CO2=True, CO=True, H2O=True, NH3=True,
                          O2=True, O3=True, C2H2=True, C2H4=True, C2H6=True,
                          H2CO=True, H2S=True, HCl=True, HCN=True, HF=True,
                          MgH=True, N2=True, NO=True, NO2=True, OCS=True,
                          OH=True, PH3=True, SH=True, SiH=True, SiO=True,
                          SO2=True, TiO=True, VO=True, Na=True, K=True,
                          scattering=True, CIA=True, **kwargs):

    with open(os.path.join(os.path.dirname(__file__),
                           'include/selectChem.in')) as f:
        selectchem_template = f.read()

    contents = selectchem_template.format(int(CH4), int(CO2), int(CO), int(H2O),
                                          int(NH3), int(O2), int(O3), int(C2H2),
                                          int(C2H4), int(C2H6), int(H2CO),
                                          int(H2S), int(HCl), int(HCN), int(HF),
                                          int(MgH), int(N2), int(NO), int(NO2),
                                          int(OCS), int(OH), int(PH3), int(SH),
                                          int(SiH), int(SiO), int(SO2), int(TiO),
                                          int(VO), int(Na), int(K),
                                          int(scattering), int(CIA))

    with open('selectChem.in', 'w') as f:
        f.write(contents)

def check_other_input(base_dir=None, opac_CH4='/Opac/opacCH4.dat', opac_C2H2='/Opac/opacC2H2.dat',
                       opac_C2H4='/Opac/opacC2H4.dat', opac_C2H6='/Opac/opacC2H6.dat',
                       opac_CO='/Opac/opacCO.dat', opac_CO2='/Opac/opacCO2.dat',
                       opac_H2CO='/Opac/opacH2CO.dat', opac_H2O='/Opac/opacH2O.dat',
                       opac_H2S='/Opac/opacH2S.dat', opac_HCN='/Opac/opacHCN.dat',
                       opac_HCl='/Opac/opacHCl.dat', opac_HF='/Opac/opacHF.dat',
                       opac_MgH='/Opac/opacMgH.dat', opac_N2='/Opac/opacN2.dat',
                       opac_NH3='/Opac/opacNH3.dat', opac_NO='/Opac/opacNO.dat',
                       opac_NO2='/Opac/opacNO2.dat', opac_O2='/Opac/opacO2.dat',
                       opac_O3='/Opac/opacO3.dat', opac_OCS='/Opac/opacOCS.dat',
                       opac_OH='/Opac/opacOH.dat', opac_PH3='/Opac/opacPH3.dat',
                       opac_SH='/Opac/opacSH.dat', opac_SO2='/Opac/opacSO2.dat',
                       opac_SiH='/Opac/opacSiH.dat', opac_SiO='/Opac/opacSiO.dat',
                       opac_TiO='/Opac/opacTiO.dat', opac_VO='/Opac/opacVO.dat',
                       opac_Na='/Opac/opacNa.dat', opac_K='/Opac/opacK.dat',
                       opac_CIA='/Opac/opacCIA.dat',
                       NT=30, Tmin=100.0, Tmax=3000.0,
                       NP=13, Pmin=1.0e-4, Pmax=1.0e8, Nlambda=4616, Ntau=334, **kwargs):

    sig = inspect.signature(create_other_input)
    for param in sig.parameters.keys():
        if 'opac' in param:
            if not os.path.isfile(base_dir + param):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                                base_dir + os.path.dirname(output_file))

    int(NT)
    int(NP)
    int(Nlambda)
    int(Ntau)
    float(Tmin)
    float(Tmax)
    float(Pmin)
    float(Pmax)

def create_other_input(opac_CH4='/Opac/opacCH4.dat', opac_C2H2='/Opac/opacC2H2.dat',
                       opac_C2H4='/Opac/opacC2H4.dat', opac_C2H6='/Opac/opacC2H6.dat',
                       opac_CO='/Opac/opacCO.dat', opac_CO2='/Opac/opacCO2.dat',
                       opac_H2CO='/Opac/opacH2CO.dat', opac_H2O='/Opac/opacH2O.dat',
                       opac_H2S='/Opac/opacH2S.dat', opac_HCN='/Opac/opacHCN.dat',
                       opac_HCl='/Opac/opacHCl.dat', opac_HF='/Opac/opacHF.dat',
                       opac_MgH='/Opac/opacMgH.dat', opac_N2='/Opac/opacN2.dat',
                       opac_NH3='/Opac/opacNH3.dat', opac_NO='/Opac/opacNO.dat',
                       opac_NO2='/Opac/opacNO2.dat', opac_O2='/Opac/opacO2.dat',
                       opac_O3='/Opac/opacO3.dat', opac_OCS='/Opac/opacOCS.dat',
                       opac_OH='/Opac/opacOH.dat', opac_PH3='/Opac/opacPH3.dat',
                       opac_SH='/Opac/opacSH.dat', opac_SO2='/Opac/opacSO2.dat',
                       opac_SiH='/Opac/opacSiH.dat', opac_SiO='/Opac/opacSiO.dat',
                       opac_TiO='/Opac/opacTiO.dat', opac_VO='/Opac/opacVO.dat',
                       opac_Na='/Opac/opacNa.dat', opac_K='/Opac/opacK.dat',
                       opac_CIA='/Opac/opacCIA.dat',
                       NT=30, Tmin=100.0, Tmax=3000.0,
                       NP=13, Pmin=1.0e-4, Pmax=1.0e8, Nlambda=4616, Ntau=334,
                       **kwargs):


    with open(os.path.join(os.path.dirname(__file__),
                           'include/otherInput.in')) as f:
        otherinput_template = f.read()

    contents = otherinput_template.format(opac_CH4, opac_C2H2, opac_C2H4,
                                          opac_C2H6, opac_CO, opac_CO2,
                                          opac_H2CO, opac_H2O, opac_H2S,
                                          opac_HCN, opac_HCl, opac_HF, opac_MgH,
                                          opac_N2, opac_NH3, opac_NO, opac_NO2,
                                          opac_O2, opac_O3, opac_OCS, opac_OH,
                                          opac_PH3, opac_SH, opac_SO2, opac_SiH,
                                          opac_SiO, opac_TiO, opac_VO, opac_Na,
                                          opac_K, opac_CIA, NT, Tmin, Tmax, NP,
                                          Pmin, Pmax, Nlambda, Ntau)

    with open('otherInput.in', 'w') as f:
        f.write(contents)

def download_url_dir(base, dirname):
    """
    Download all of the contents of a directory from GitHub.
    """
    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass

    r = requests.get(urljoin(base, dirname))
    for entry in r.json():
        print("Downloading {} file {}".format(dirname, entry['name']))
        text = requests.get(entry['download_url']).text
        with open(os.path.join(dirname, entry['name']), 'w') as f:
            f.write(text)

def collect_exotransmit_data():
    """
    Get all of the Exo_Transmit ancillary data from the original Github repo.
    """
    base_url = 'https://api.github.com/repos/ExoCTK/exoctk_data/contents/exotransmit/'
    download_url_dir(base_url, 'EOS')
    download_url_dir(base_url, 'Opac')
    download_url_dir(base_url, 'T_P')
    r = requests.get(base_url)
    for entry in r.json():
        if entry['name'].endswith('.in'):
            print("Downloading {} file".format(entry['name']))
            text = requests.get(entry['download_url']).text
            if entry['name'] == 'userInput.in':
                print('Replacing /YOUR_PATH with {} in userInput.in'.format(os.path.abspath(os.curdir)))
                text = text.replace('/YOUR_PATH', os.path.abspath(os.curdir))
            with open(entry['name'], 'w') as f:
                f.write(text)
    try:
        os.mkdir('Spectra')
    except FileExistsError:
        pass