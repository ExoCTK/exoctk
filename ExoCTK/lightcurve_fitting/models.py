"""Base and child classes to handle models
used to fit light curves

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as q
import batman
import copy

from .parameters import Parameters


class Model:
    def __init__(self, **kwargs):
        """
        Create a model instance
        """
        # Set up model attributes
        self.name = 'New Model'
        self._time = None
        self._flux = None
        self._units = q.day
        self._parameters = None
        self.components = None
        
        # Store the arguments as attributes
        for arg, val in kwargs.items():
            setattr(self, arg, val)


    def __mul__(self, other):
        """Multiply model components to make a combined model
        
        Parameters
        ----------
        other: ExoCTK.lightcurve_fitting.models.Model
            The model to multiply
        
        Returns
        -------
        ExoCTK.lightcurve_fitting.lightcurve.Model
            The combined model
        """
        # Make sure it is the right type
        if not all([hasattr(other, attr) for attr in ['units','flux','time']]):
            raise TypeError('Only another Model instance may be multiplied.')
        
        # Convert other time axis to same units
        other.units = self.units
        
        return CompositeModel([copy.copy(self), other])
        
        
    def interp(self, new_time, units=q.day):
        """Interpolate the flux to a new time axis
        
        Parameters
        ----------
        new_time: sequence, astropy.units.quantity.Quantity
            The time array
        units: str, astropy.units.core.Unit, astropy.units.core.IrreducibleUnit
            The units of the input time_array, 'day' by default
        """
        # Check the type
        if not isinstance(new_time, (np.ndarray, tuple, list, q.quantity.Quantity)):
            raise TypeError("Time axis must be a tuple, list, astropy quantity, or numpy array.")
        
        # Use given units if provided
        if hasattr(new_time, 'unit'):
            units = new_time.unit
            new_time = new_time.value
            
        # Calculate the new_time
        new_time = (np.array(new_time)*units).to(self.units).value
        
        # Calculate the new flux
        self.flux = np.interp(new_time, self.time, self.flux)
        
        # Set the new time axis
        self.time = new_time
        
        
    @property
    def parameters(self):
        """A getter for the parameters"""
        return self._parameters


    @parameters.setter
    def parameters(self, params):
        """A setter for the parameters"""
        # Process if it is a parameters file
        if isinstance(params, str) and os.file.exists(params):
            params = Parameters(params)
            
        # Or a Parameters instance
        if not isinstance(params, (Parameters, type(None))):
            raise TypeError("'params' argument must be a JSON file, ascii file, or ExoCTK.lightcurve_fitting.parameters.Parameters instance.")
            
        # Set the parameters attribute
        self._parameters = params
        
        
    @property
    def time(self):
        """A getter for the time"""
        return self._time
        
        
    @time.setter
    def time(self, time_array, units=q.day):
        """A setter for the time
        
        Parameters
        ----------
        time_array: sequence, astropy.units.quantity.Quantity
            The time array
        units: str, astropy.units.core.Unit, astropy.units.core.IrreducibleUnit
            The units of the input time_array, 'day' by default
        """
        # Check the type
        if not isinstance(time_array, (np.ndarray, tuple, list, q.quantity.Quantity)):
            raise TypeError("Time axis must be a tuple, list, astropy quantity, or numpy array.")
        
        # Use given units if provided
        if hasattr(time_array, 'unit'):
            units = time_array.unit
            time_array = time_array.value
        
        # Set the array
        self._time = (np.array(time_array)*units).to(self.units).value
        
        
    @property
    def units(self):
        """A getter for the units"""
        return self._units
        
        
    @units.setter
    def units(self, units):
        """A setter for the units
        
        Parameters
        ----------
        units: str, astropy.units.core.Unit, astropy.units.core.IrreducibleUnit
            The time units
        """
        # Convert string to unit
        if isinstance(units, str):
            units = q.Unit(units)
            
        # Check the type
        if not isinstance(units, (q.core.IrreducibleUnit, q.core.Unit)):
            raise TypeError("units axis must be a tuple, list, astropy quantity, or numpy array.")
            
        # Make sure they are time units
        _ = units.to(q.day)
            
        # Set the attribute
        self._units = units
        
        # Update the time
        if self.time is not None:
            scale = self.units.to(units)
            self.time = self.time*scale
        
        
    @property
    def flux(self):
        """A getter for the flux"""
        return self._flux
        
        
    @flux.setter
    def flux(self, flux_array):
        """A setter for the flux
        
        Parameters
        ----------
        flux_array: sequence
            The flux array
        """
        # Check the type
        if not isinstance(flux_array, (np.ndarray, tuple, list)):
            raise TypeError("flux axis must be a tuple, list, or numpy array.")
        
        # Set the array
        self._flux = np.array(flux_array)
        
        
    def plot(self, time, components=False, **kwargs):
        """Plot the model
        
        Parameters
        ----------
        components: bool
            Plot all model components
        """
        # Set the time
        self.time = time
        
        flux = self.eval(**kwargs)
        plt.plot(self.time, flux, label=self.name)
        
        if components and self.components is not None:
            for comp in self.components:
                flux = comp.eval(**kwargs)
                plt.plot(self.time, flux, label=comp.name)
            
        plt.xlabel(self.units.long_names[0])
        plt.ylabel('Flux')
        
        plt.legend(loc=0)


class CompositeModel(Model):
    """A class to create composite models"""
    def __init__(self, models, **kwargs):
        """Initialize the composite model
        
        Parameters
        ----------
        models: sequence
            The list of models
        """
        # Inherit from Model calss
        super().__init__(**kwargs)
        
        # Store the models
        self.components = models
        
    def eval(self, **kwargs):
        """Evaluate the model components"""
        # Get the time
        if self.time is None:
            self.time = kwargs.get('time')
        
        # Empty flux
        flux = 1.
        
        # Evaluate flux at each model
        for model in self.components:
            flux *= model.eval(**kwargs)
            
        return flux


class PolynomialModel(Model):
    """Polynomial Model"""
    def __init__(self, **kwargs):
        """Initialize the polynomial model
        """
        # Inherit from Model calss
        super().__init__(**kwargs)
        
        # Check for Parameters instance
        self.parameters = kwargs.get('parameters')
            
        # Generate parameters from kwargs if necessary
        if self.parameters is None:
            coeffs = self._parse_coeffs(kwargs)
            params = {'c{}'.format(n): coeff for n,coeff in enumerate(coeffs[::-1])}
            self.parameters = Parameters(**params)
            
            
    @staticmethod
    def _parse_coeffs(coeff_dict):
        """Convert dict of 'c#' coefficients into a list
        of coefficients in decreasing order, i.e. ['c2','c1','c0']
        
        Parameters
        ----------
        coeff_dict: dict
            The dictionary of coefficients
        
        Returns
        -------
        np.ndarray
            The sequence of coefficient values
        """
        # Parse 'c#' keyword arguments as coefficients
        coeffs = np.zeros(10)
        for k,v in coeff_dict.items():
            if k.lower().startswith('c') and k[1:].isdigit():
                coeffs[int(k[1:])] = v
           
        # Trim zeros and reverse
        return np.trim_zeros(coeffs)[::-1]
        
        
    def eval(self, **kwargs):
        """Evaluate the function with the given values"""
        # Get the time
        if self.time is None:
            self.time = kwargs.get('time')
        
        coeffs = [coeff[1] for coeff in self.parameters.list][::-1]
        
        # Create the polynomial from the coeffs
        poly = np.poly1d(coeffs)
        
        # Evaluate the polynomial
        return np.polyval(poly, self.time)


class TransitModel(Model):
    """Transit Model"""
    def __init__(self, **kwargs):
        """Initialize the transit model
        """
        # Inherit from Model calss
        super().__init__(**kwargs)
        
        # Check for Parameters instance
        self.parameters = kwargs.get('parameters')
        
        # Generate parameters from kwargs if necessary
        if self.parameters is None:
            self.parameters = Parameters(**kwargs)
        
        
    def eval(self, **kwargs):
        """Evaluate the function with the given values"""
        # Get the time
        if self.time is None:
            self.time = kwargs.get('time')
        
        # Generate with batman
        bm_params = batman.TransitParams()
        
        # Set all parameters
        for p in self.parameters.list:
            setattr(bm_params, p[0], p[1])
            
        # Set t0 without units
        bm_params.t0 = self.parameters.t0.value*q.day.to(self.units)
        
        # Combine limb darkening coeffs
        bm_params.u = [self.parameters.u1.value, self.parameters.u2.value]
        
        # Make the eclipse
        m_eclipse = batman.TransitModel(bm_params, self.time, transittype=self.parameters.transittype.value)

        # OoT == Out of transit    
        OoT_curvature = self.parameters.offset.value+self.parameters.slope.value*(self.time-self.time.mean())+self.parameters.curvature.value*(self.time-self.time.mean())**2

        # Evaluate the light curve
        return m_eclipse.light_curve(bm_params) * OoT_curvature

