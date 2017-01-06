from distutils.extension import Extension

def get_extensions():
    return [Extension(name='ExoCTK.pal.exotransmit', sources=['ExoCTK/pal/exotransmit.pyx'],
                     include_dirs=['numpy', 'ExoCTK/pal/include'])]

def get_package_data():
    return {'ExoCTK.pal': ['data/Opac/*', 'data/EOS/*', 'data/T_P/*', 'include/*']}
