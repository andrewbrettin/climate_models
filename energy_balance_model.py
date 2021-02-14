import numpy as np
import xarray as xr

from scipy.constants import sigma
from scipy.integrate import solve_ivp


class EnergyBalanceModel:
    """
    Energy balance model is meant to solve the 1D differential equation
    
    .. math::
        c_p \frac{\partial T(\varphi)}{\partial t} = S_0(\alpha) (1 - \alpha(\phi))
        - \left(1 - \frac{\epsilon}{2}\right) \sigma T(\varphi)^4 + 
        + D \frac{\partial^2 T(\varphi)}{\partial \varphi^2}
    
    
    Parameters: 
        c_p : float
            Heat capacity per square meter [j K^-1 m^-2].
        S_0 : callable
            Insolation per square meter (function of latitude) [W m^-2].
        alpha : callable
            Albedo as a function of the latitude phi [unitless].
        epsilon : float
            Emissivity of atmosphere [unitless].
        D : float
            Diffusivity [W].
        grid : xarray.DataArray
            Temperature grid as a function of time and latitude. Initialized as an empty array,
            and populated using EnergyBalanceModel.run_model().
    Methods:
        run_model
    """
    
    
    def __init__(self):
        # To do
        self.c_p = 1
        self.S_0 = lambda phi : 1 # To do
        self.alpha = lambda phi : 0.6 if abs(phi) < 273 else 0.3
        self.epsilon = 0.76
        self.D = 1
        self.grid = xr.DataArray([])
        
    def __init__(self, c_p=4e6, S_0=lambda phi : 1,
                 alpha=lambda phi: 0.6 if abs(phi) < 273 else 0.3,
                 epsilon=0.76, D=6e15):
        self.c_p = c_p
        self.S_0 = S_0
        self.alpha = alpha
        self.epsilon = epsilon
        self.D = D
        self.grid = xr.DataArray([])
    
    def run_model(self, res=1.0, dt=100, N=1000, 
                  T_0=lambda x : 300 - (5/810) * x**2, method='RK45'):
        """
        Timesteps an EnergyBalanceModel. Uses xarray to handle the grid and scipy.integrate to timetep.
        
        Parameters:
            res : float
                Degree resolution [degrees]. Default is 1.0; accepts values 
                0.1, 0.5, 1, 2, 5, 10, 15, 30, and 45.
            dt : float
                Timestep of simulation [years]
            N : int
                Number of timesteps to integrate.
            T_0 : callable
                Function of latitude giving the initial temperature profile [K].
            method : string
                Timestep method. Accepts 'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', or 'LSODA'
            
        """
        
        if not (res in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 30.0, 45.0]):
            raise ValueError('Unacceptable degree resolution. \n \
                             Acceptable values include 0.1, 0.5, 1, 2, 5, 10, 15, 30, and 45.')
        if not (method in ['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA']):
            raise ValueError('ODE timestepping method not defined.')
            
        # Initialize grid
        ts = dt*np.arange(0,N)
        phis = np.arange(-90, 90+res, res)
        
        data = np.zeros((N, len(phis)))
        self.grid = xr.DataArray(data, dims=('time', 'latitude'), coords={'time': ts, 'latitude': phis})
        
        grid.loc[{'time': 0}] = [T_0(phi) for phi in phis]
        
        # Define RHS of 
        
        # Timestep grid
        for i in range(1, N):
            # RK45 timestep:
            solve_ivp
            
        
        
        return T_grid
        
        
        
        
        
        
        
        
        