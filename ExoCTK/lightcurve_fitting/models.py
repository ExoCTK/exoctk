"""Base and child classes to handle models
used to fit light curves

Author: Joe Filippazzo
Email: jfilippazzo@stsci.edu
"""
import numpy as np
import matplotlib.pyplot as plt

class Model:
    def __init__(self, model_name, time=None, position=None, **kwargs):
        """
        Initialize a supported model

        Parameters
        ----------
        model_name: str
            The name of the model
        model_dir: 
        """
        # Store the name
        self.name = model_name
        
        self.model = None


    def __multiply__(self, other):
        """Multiply model components to make a combined model
        
        Parameters
        ----------
        other: ExoCTK.tlc.lightcurve.Model
            The model to multiply
        
        Returns
        -------
        ExoCTK.tlc.lightcurve.Model
            The combined model
        """
        # Make sure it is the right type
        if not isinstance(other, Model):
            raise TypeError('Only another Model instance may be multiplied.')

        new_model = self.model*other.model

        return CustomModel(new_model)
        
    
    def plot(self):
        """Plot the model"""
        plt.figure()
        
        if self.model is not None:
            plt.plot(self.model)


class CustomModel(Model):
    """Custom Model class, which can be passed arbitrary arrays"""
    def __init__(self, model_name, model):
        """Initialize the custom model"""
        super().__init__(model_name)

        # Set the model
        if not isinstance(model, np.ndarray):
            raise TypeError('Bad model: must be a sequence.')

        self.model = np.array(model)


class LinearModel(Model):
    """Linear Model"""
    def __init__(self, model_name, value=1):
        """Initialize the linear model"""
        super().__init__(model_name)

        self.model = np.ones(100)*value


class QuadraticModel(Model):
    """Quadratic Model"""
    def __init__(self, model_name):
        """Initialize the quadratic model"""
        super().__init__(model_name)
        
        self.model = np.arange(100)**2


class StellarModel(Model):
    """Stellar Model"""
    def __init__(self, model_name):
        """Initialize the stellar model"""
        super().__init__(model_name)

        self.model = np.arange(100)**2


class SystematicModel(Model):
    """Systematic Model"""
    def __init__(self, model_name):
        """Initialize the systematic model"""
        super().__init__(model_name)

        self.model = np.arange(100)**2


class TransitModel(Model):
    """Transit Model"""
    def __init__(self, model_name):
        """Initialize the transit model"""
        super().__init__(model_name)

        self.model = np.arange(100)**2

