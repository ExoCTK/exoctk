from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(ext_modules=[Extension("include/_tran_module", ["include/_tau.c", "include/_tau_wrap.c", "include/_xsects.c"])], 
      include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs(),)