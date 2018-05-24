from distutils.extension import Extension

def get_extensions():
    return [Extension(name='ExoCTK.forward_models._chimera._tran_module',
        sources=["ExoCTK/forward_models/_chimera/src/_tau.c",
            "ExoCTK/forward_models/_chimera/src/_tau_wrap.c",
            "ExoCTK/forward_models/_chimera/src/_xsects.c"],
                     include_dirs=['numpy', 'ExoCTK/forward_models/_chimera/include'])]

def get_package_data():
    return {'ExoCTK.forward_models._chimera': ['include/*', 'data/*']}
