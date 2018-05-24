from distutils.extension import Extension
import os

def get_extensions():
    cfiles = []
    for f in os.listdir('ExoCTK/forward_models/include/'):
        if f.endswith('.c') and f != 'main_transmission.c':
            cfiles.append(os.path.join('ExoCTK/forward_models/include/', f))

    return [Extension(name='ExoCTK.pal._exotransmit_wrapper',
                      sources=['ExoCTK/forward_models/_exotransmit_wrapper.pyx']+cfiles,
                      include_dirs=['numpy', 'ExoCTK/forward_models/include'])]

def get_package_data():
    return {'ExoCTK.forward_models': ['include/*']}
