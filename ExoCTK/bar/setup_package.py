# from distutils.extension import Extension
# 
# def get_extensions():
#     return [Extension(name='ExoCTK.pal.exotransmit', sources=['ExoCTK/pal/exotransmit.pyx'],
#                      include_dirs=['numpy', 'ExoCTK/pal/include'])]
# 
# def get_package_data():
#     return {'ExoCTK.pal': ['data/Opac/*', 'data/EOS/*', 'data/T_P/*', 'include/*']}
# 

from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(ext_modules=[Extension("_tran_module", ["_tau.c", "_tau_wrap.c", "_xsects.c"])], 
      include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs(),)
