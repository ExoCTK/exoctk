from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(ext_modules=[Extension("include/_tran_module", ["_tau.c", "_tau_wrap.c", "_xsects.c"])], 
      include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs(),)