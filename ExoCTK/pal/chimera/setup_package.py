from distutils.extension import Extension

def get_extensions():
    return [Extension(name='ExoCTK.pal.chimera._tran_module', sources=["ExoCTK/pal/chimera/src/_tau.c", "ExoCTK/pal/chimera/src/_tau_wrap.c", "ExoCTK/pal/chimera/src/_xsects.c"],
                     include_dirs=['numpy', 'ExoCTK/pal/chimera/include'])]

def get_package_data():
    return {'ExoCTK.pal.chimera': ['include/*']}
