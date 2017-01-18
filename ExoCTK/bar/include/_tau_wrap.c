//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include "_tau.h"
#include "_xsects.h"

static char module_docstring[] = "Forward model.";

static char _calc_tau_docstring[] = "Returns tau at each level and wavenumber.";
static char _init_xsects_docstring[] ="Reads in precomputed cross-sections and rebins them according to input T and P grids";

static PyObject * _tau_wrap(PyObject* self, PyObject* args);
static PyObject * _init_xsects_wrap(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"_tau_wrap", _tau_wrap, METH_VARARGS, _calc_tau_docstring},
  {"_init_xsects_wrap", _init_xsects_wrap, METH_VARARGS, _init_xsects_docstring},
  {NULL}
};

#if PY_MAJOR_VERSION >= 3
// C API changes for Python 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_tran_module",     /* m_name */
    module_docstring,  /* m_doc */
    -1,                  /* m_size */
    module_methods,    /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};
PyMODINIT_FUNC PyInit__tran_module(void) {
  PyObject* m = PyModule_Create(&moduledef);
  if (m == NULL) 
    return NULL;
  import_array();
  return m;
}
#else
PyMODINIT_FUNC init_tran_module(void) {
  PyObject* m = Py_InitModule3("_tran_module", module_methods,
			       module_docstring);
  if (m == NULL) return;

  import_array();
}
#endif
/*
Function: _tau_wrap
-------------
Arguments: Xsects, RXsects, Z, Pavg, Tavg, Fractions, Tau, r0
*/
static PyObject* _tau_wrap(PyObject* self, PyObject* args) {
  double r0;
  PyObject *Xsect_obj, *RXsect_obj, *Z_obj, *P_obj, *T_obj, *F_obj, *Tau_obj;

  if (!PyArg_ParseTuple(args, "OOOOOOOd", &Xsect_obj,
			&RXsect_obj, &Z_obj, &P_obj,
			&T_obj, &F_obj, &Tau_obj, &r0)) {
      return NULL;
  }

  PyObject* Xsect_arr = PyArray_FROM_OTF(Xsect_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* RXsect_arr = PyArray_FROM_OTF(RXsect_obj,
					  NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Z_arr = PyArray_FROM_OTF(Z_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* P_arr = PyArray_FROM_OTF(P_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* T_arr = PyArray_FROM_OTF(T_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* F_arr = PyArray_FROM_OTF(F_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Tau_arr = PyArray_FROM_OTF(Tau_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);

  //If any of the arrays are unable to be created,
  //PyArray_FROM_OTF will set the error status and 
  //we return NULL indicating there was an error.
  if (!(Xsect_arr && RXsect_arr && Z_arr
	&& P_arr && T_arr && F_arr && Tau_arr)) {
    Py_XDECREF(Xsect_arr);
    Py_XDECREF(RXsect_arr);
    Py_XDECREF(Z_arr);
    Py_XDECREF(P_arr);
    Py_XDECREF(T_arr);
    Py_XDECREF(F_arr);
    Py_XDECREF(Tau_arr);
    return NULL;
  }

  int numlevels = (int) PyArray_DIM(Xsect_arr, 0) + 1;
  int numbins = (int) PyArray_DIM(Xsect_arr, 1);
  int numsources = (int) PyArray_DIM(Xsect_arr, 2);

  double kb = 1.38e-23;
  
  double *Xsects = (double *)PyArray_DATA(Xsect_arr);
  double *RXsects = (double *)PyArray_DATA(RXsect_arr);
  double *Z = (double *)PyArray_DATA(Z_arr);
  double *Pavg = (double *)PyArray_DATA(P_arr);
  double *Tavg = (double *)PyArray_DATA(T_arr);
  double *Fractions = (double *)PyArray_DATA(F_arr);
  double *Tau = (double *)PyArray_DATA(Tau_arr);

  //printf("The shape of tau is: (%d, %d)\n", numbins, numlevels);
  //printf("The first values of tau are: %f, %f, %f\n", Tau[0], Tau[1], Tau[2]);
  _calc_tau(Xsects, RXsects, Z, Pavg, Tavg, Fractions, numlevels, numbins, numsources, Tau, kb, r0);

  //Dont forget to clean up by decreasing reference count on 
  //arrays so their memory can be deallocated
  Py_DECREF(Xsect_arr);
  Py_DECREF(RXsect_arr);
  Py_DECREF(Z_arr);
  Py_DECREF(P_arr);
  Py_DECREF(T_arr);
  Py_DECREF(F_arr);
  Py_DECREF(Tau_arr);
  
  //Return None
  PyObject *ret = Py_BuildValue("");
  return ret;
}

static PyObject* _init_xsects_wrap(PyObject* self, PyObject* args) {
  int bin_offset;
  double scatter_power, scatter_coeff;
  PyObject *Xsectout_obj, *RXsectout_obj, *Xsectin_obj, *Tgrid_obj, 
    *Pgrid_obj, *Tavg_obj, *Pavg_obj, *wno_obj;

  if (!PyArg_ParseTuple(args, "OOOOOOOOidd", &Xsectout_obj,
			&RXsectout_obj, &Xsectin_obj, &Tgrid_obj,
			&Pgrid_obj, &Tavg_obj, &Pavg_obj,
			&wno_obj, &bin_offset, &scatter_coeff, &scatter_power)) {
      return NULL;
  }

  PyObject* Xsectout_arr = PyArray_FROM_OTF(Xsectout_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* RXsectout_arr = PyArray_FROM_OTF(RXsectout_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Xsectin_arr = PyArray_FROM_OTF(Xsectin_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Tgrid_arr = PyArray_FROM_OTF(Tgrid_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Pgrid_arr = PyArray_FROM_OTF(Pgrid_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Tavg_arr = PyArray_FROM_OTF(Tavg_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* Pavg_arr = PyArray_FROM_OTF(Pavg_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);
  PyObject* wno_arr = PyArray_FROM_OTF(wno_obj, NPY_DOUBLE,
				     NPY_IN_ARRAY);

  //If any of the arrays are unable to be created,
  //PyArray_FROM_OTF will set the error status and 
  //we return NULL indicating there was an error.
  if (!(Xsectout_arr && RXsectout_arr && Xsectin_arr
	&& Tgrid_arr && Pgrid_arr && Tavg_arr && Pavg_arr
	&& wno_arr)) {
    Py_XDECREF(Xsectout_arr);
    Py_XDECREF(RXsectout_arr);
    Py_XDECREF(Xsectin_arr);
    Py_XDECREF(Tgrid_arr);
    Py_XDECREF(Pgrid_arr);
    Py_XDECREF(Tavg_arr);
    Py_XDECREF(Pavg_arr);
    Py_XDECREF(wno_arr);
    return NULL;
  }

  int numgases = (int) PyArray_DIM(Xsectout_arr, 2);
  int numlevels = (int) PyArray_DIM(Xsectout_arr, 0) + 1;
  int numbins = (int) PyArray_DIM(Xsectout_arr, 1);
  int numprebins = (int) PyArray_DIM(Xsectin_arr, 3);
  int Tgridlen = (int) PyArray_DIM(Tgrid_arr, 0);
  int Pgridlen = (int) PyArray_DIM(Pgrid_arr, 0);

  double *Xsectsout = (double *)PyArray_DATA(Xsectout_arr);
  double *RXsectsout = (double *)PyArray_DATA(RXsectout_arr);
  double *Xsectsin = (double *)PyArray_DATA(Xsectin_arr);
  double *Tgrid = (double *)PyArray_DATA(Tgrid_arr);
  double *Pgrid = (double *)PyArray_DATA(Pgrid_arr);
  double *Tavg = (double *)PyArray_DATA(Tavg_arr);
  double *Pavg = (double *)PyArray_DATA(Pavg_arr);
  double *wno = (double *)PyArray_DATA(wno_arr);

  _init_xsects(Xsectsout, RXsectsout, Xsectsin, Tgrid, Pgrid, Tavg,
	       Pavg, wno, scatter_coeff, scatter_power, numlevels,
	       numbins, bin_offset, numprebins, numgases, Tgridlen, Pgridlen);

  Py_DECREF(Xsectout_arr);
  Py_DECREF(RXsectout_arr);
  Py_DECREF(Xsectin_arr);
  Py_DECREF(Tgrid_arr);
  Py_DECREF(Pgrid_arr);
  Py_DECREF(Tavg_arr);
  Py_DECREF(Pavg_arr);
  Py_DECREF(wno_arr);

  PyObject *ret = Py_BuildValue("");
  return ret;
}
