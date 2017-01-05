from distutils.extension import Extension

def get_extensions():
    return [Extension(name='ExoCTK.exotransmit', sources=['ExoCTK/exotransmit/exotransmit.pyx'],
                     include_dirs=['numpy', 'ExoCTK/exotransmit/include'])]