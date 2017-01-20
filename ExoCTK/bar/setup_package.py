from distutils.extension import Extension

def get_extensions():
    return [Extension(name='ExoCTK.bar._tran_module', sources=["ExoCTK/bar/include/_tau.c", "ExoCTK/bar/include/_tau_wrap.c", "ExoCTK/bar/include/_xsects.c"],
                     include_dirs=['numpy', 'ExoCTK/bar/include'])]

def get_package_data():
    return {'ExoCTK.bar': ['data/*', 'abscoeff/*', 'chains/*', 'include/*']}
