#! /usr/bin/env python

"""Tests for the ``lightcurve_fitting`` package.

Authors
-------

    Joe Filippazzo

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_lightcurve_fitting.py
"""

import unittest

import numpy as np

from ..lightcurve_fitting import lightcurve, models, parameters, simulations


class TestLightcurve(unittest.TestCase):
    """Tests for the lightcurve.py module"""
    def setUp(self):
        """Setup for the lightcurve"""
        self.time = np.linspace(0, 1, 100)
        self.unc = np.random.uniform(low=1E-4, high=0.01, size=100)
        self.flux = np.random.normal(np.ones_like(self.time), scale=self.unc)

    def test_lightcurve(self):
        """Test that a LightCurve object can be created"""
        self.lc = lightcurve.LightCurve(self.time, self.flux, self.unc, name='Data')

        # Test that parameters can be assigned
        lin1 = models.PolynomialModel(c1=0.0005, c0=0.997, name='linear 1')
        lin2 = models.PolynomialModel(c1=0.001, c0=0.92, name='linear 2')
        comp_model = lin1*lin2

        # Test the fitting routine
        self.lc.fit(comp_model)


class TestModels(unittest.TestCase):
    """Tests for the models.py module"""
    def setUp(self):
        """Setup for the tests"""
        # Set time to use for evaluations
        self.time = np.linspace(0, 1, 100)

    def test_model(self):
        """Tests for the generic Model class"""
        # Test model creation
        name = 'Model 1'
        self.model = models.Model(name=name)
        self.assertEqual(self.model.name, name)

        # Test model units
        self.assertEqual(str(self.model.units), 'd')
        self.model.units = 'MJD'
        self.assertEqual(self.model.units, 'MJD')
        self.assertRaises(TypeError, setattr, self.model.units, 'foobar')

    def test_compositemodel(self):
        """Tests for the CompositeModel class"""
        model1 = models.Model()
        model2 = models.Model()
        self.comp_model = model1*model2
        self.comp_model.name = 'composite'

    def test_polynomialmodel(self):
        """Tests for the PolynomialModel class"""
        # Create the model
        self.lin_model = models.PolynomialModel(c1=0.0005, c0=0.997, name='linear')

        # Evaluate and test output
        self.lin_model.time = self.time
        vals = self.lin_model.eval()
        self.assertEqual(vals.size, self.time.size)

    def test_transitmodel(self):
        """Tests for the TransitModel class"""
        # Set the intial parameters
        params = parameters.Parameters()
        params.rp = 0.22, 'free', 0.0, 0.4  # rprs
        params.per = 10.721490, 'fixed'
        params.t0 = 0.48, 'free', 0, 1
        params.inc = 89.7, 'free', 80., 90.
        params.a = 18.2, 'free', 15., 20.    # aprs
        params.ecc = 0., 'fixed'
        params.w = 90., 'fixed'             # omega
        params.limb_dark = '4-parameter', 'independent'
        params.transittype = 'primary', 'independent'
        params.u1 = 0.1, 'free', 0., 1.
        params.u2 = 0.1, 'free', 0., 1.
        params.u3 = 0.1, 'free', 0., 1.
        params.u4 = 0.1, 'free', 0., 1.

        # Make the transit model
        self.t_model = models.TransitModel(parameters=params, name='transit')


class TestParameters(unittest.TestCase):
    """Tests for the parameters.py module"""
    def setUp(self):
        """Setup for the tests"""
        pass

    def test_parameter(self):
        """Test that a Parameter object can be created"""
        # Create the parameter
        pname = 'p1'
        pval = 12.34
        ptype = 'free'
        pmn = 10
        pmx = 15
        self.param = parameters.Parameter(pname, pval, ptype, pmn, pmx)

        # Test bogus input
        self.assertRaises(TypeError, parameters.Parameter, 123)
        self.assertRaises(ValueError, parameters.Parameter, 'foo', 123, 123)

        # Test the attributes
        self.assertEqual(self.param.name, pname)
        self.assertEqual(self.param.value, pval)
        self.assertEqual(self.param.ptype, ptype)
        self.assertEqual(self.param.mn, pmn)
        self.assertEqual(self.param.mx, pmx)
        self.assertEqual(self.param.values, (pname, pval, ptype, pmn, pmx))

    def test_parameters(self):
        """Test that a Parameters object can be created"""
        self.params = parameters.Parameters()
        self.params.param1 = 123.456
        self.params.param2 = 234.567, True, 200, 300

        # Test the auto attribute assignment
        self.assertEqual(self.params.param1.values, ('param1', 123.456, 'free'))
        self.assertEqual(self.params.param2.values, ('param2', 234.567, 'free', 200, 300))


# class TestSimulations(unittest.TestCase):
#     """Test the simulations.py module"""
#     def setUp(self):
#         """Setup for the tests"""
#         pass
#
#     def test_simulation(self):
#         """Test the simulations can be made properly"""
#         # Test to pass
#         npts = 1234
#         time, flux, unc, params = simulations.simulate_lightcurve('WASP-19b', 0.1, npts=npts, plot=True)
#         self.assertEqual(len(time), npts)
#
#         # Test to fail
#         self.assertRaises(ValueError, simulations.simulate_lightcurve, 'foobar', 0.1)
