"""Base and child classes to handle models
used to fit light curves

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as q

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
        
        # Interpolate the other flux to the current time axis and multiply
        other.interp(self.time, self.units)
        
        # Multiply fluxes
        new_flux = self.flux*other.flux
        
        # Store the components
        components = [self, other]
        
        return Model(time=self.time, flux=new_flux, components=components)
        
        
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
        if not isinstance(params, Parameters):
            raise TypeError("'params' argument must be a JSON file, ascii file, or ExoCTK.lightcurve_fitting.lightcurve.Parameters instance.")
            
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
        
        
    def plot(self, components=False):
        """Plot the model
        
        Parameters
        ----------
        components: bool
            Plot all model components
        """
        if self.time is not None and self.flux is not None:
            plt.plot(self.time, self.flux, label=self.name)
            
            if components and self.components is not None:
                for comp in self.components:
                    comp.plot()
                
            plt.xlabel(self.units.long_names[0])
            plt.ylabel('Flux')
            
            plt.legend(loc=0)
        
        else:
            print('No model data to plot.')


class PolynomialModel(Model):
    """Polynomial Model"""
    def __init__(self, time, coeffs, **kwargs):
        """Initialize the polynomial model
        
        Parameters
        ----------
        time_array: sequence, astropy.units.quantity.Quantity
            The time array
        coeffs: float, int, sequence
            The list of coefficients, highest order first
        """
        super().__init__(time=time, **kwargs)
        
        # Create the polynomial from the coeffs
        poly = np.poly1d(coeffs)
        
        # Evaluate
        self.flux = np.polyval(poly, time)


# class QuadraticModel(Model):
#     """Quadratic Model"""
#     def __init__(self):
#         """Initialize the quadratic model"""
#         super().__init__()
#
#         self.model = np.arange(100)**2
#
#
# class StellarModel(Model):
#     """Stellar Model"""
#     def __init__(self):
#         """Initialize the stellar model"""
#         super().__init__()
#
#         self.model = np.arange(100)**2
#
#
# class SystematicModel(Model):
#     """Systematic Model"""
#     def __init__(self):
#         """Initialize the systematic model"""
#         super().__init__()
#
#         self.model = np.arange(100)**2
#
#
# class TransitModel(Model):
#     """Transit Model"""
#     def __init__(self):
#         """Initialize the transit model"""
#         super().__init__()
#
#         self.model = np.arange(100)**2

