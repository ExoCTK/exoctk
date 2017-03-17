from distutils.extension import Extension

def get_extensions():
    return [Extension(name='ExoCTK.pal._chimera._tran_module', sources=["ExoCTK/pal/_chimera/src/_tau.c", "ExoCTK/pal/_chimera/src/_tau_wrap.c", "ExoCTK/pal/_chimera/src/_xsects.c"],
                     include_dirs=['numpy', 'ExoCTK/pal/_chimera/include'])]

def get_package_data():
    return {'ExoCTK.pal._chimera': ['include/*', 'data/*']}
