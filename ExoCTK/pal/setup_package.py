from distutils.extension import Extension
import os

def get_extensions():
    cfiles = []
    for f in os.listdir('ExoCTK/pal/include/'):
        if f.endswith('.c') and f != 'main_transmission.c':
            cfiles.append(os.path.join('ExoCTK/pal/include/', f))

    return [Extension(name='ExoCTK.pal._exotransmit_wrapper',
                      sources=['ExoCTK/pal/_exotransmit_wrapper.pyx']+cfiles,
                      include_dirs=['numpy', 'ExoCTK/pal/include'])]

def get_package_data():
    return {'ExoCTK.pal': ['include/*']}
