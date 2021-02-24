"""
energy_balance_model.py

@author: andrewbrettin
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from scipy.constants import sigma
from scipy.integrate import solve_ivp
from scipy.special import legendre

import warnings
import time

class EnergyBalanceModel:
    """
    Energy balance model is meant to solve the 1D differential equation
    
    .. math::
        c_p \frac{\partial T(\varphi)}{\partial t} = S_0(\alpha) (1 - \alpha(\phi))
        - \left(1 - \frac{\epsilon}{2}\right) \sigma T(\varphi)^4 + 
        + D \frac{\partial^2 T(\varphi)}{\partial \varphi^2}
    
    
    Parameters: 
        c_p : float
            Heat capacity per square meter [J K^-1 m^-2].
        S_0 : callable
            Insolation per square meter (function of latitude) [W m^-2].
        alpha : callable
            Albedo as a function of the latitude phi [unitless].
        epsilon : float
            Emissivity of atmosphere [unitless].
        D : float
            Diffusivity [W].
        T_grid : xarray.DataArray
            Temperature grid as a function of time and latitude. Initialized as an empty array,
            and populated using EnergyBalanceModel.run_model().
    """
    
    
    def __init__(self, c_p=4e7, epsilon=0.76, D=5.30e13,
                 S_0=lambda phi : 230 * np.exp(-phi**2/1800) + 170,
                 alpha=lambda T: 0.6 if abs(T) < 273 else 0.3):
        self.c_p = c_p
        self.S_0 = np.vectorize(S_0) # Vectorized so scalar function can operate on array of phis
        self.alpha = np.vectorize(alpha)
        self.epsilon = epsilon
        self.D = D
        self.T_grid = xr.DataArray([])
    
    def run_model(self, res=1.0, dt=1e4, N=100, 
                  T_0=lambda x : 300 - (5/810) * x**2,
                 method='RK45'):
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
                Timestepping method for nondiffusive ODE component of model.
                Acceptable values include: ['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'].
            
        """
        
        if not (res in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 30.0, 45.0]):
            raise ValueError('Unacceptable degree resolution. \n \
                             Acceptable values include 0.1, 0.5, 1, 2, 5, 10, 15, 30, and 45.')
        
        if not method in ['RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA']:
            raise ValueError('ODE integration scheme undefined. See \
                scipy.integrate.solve_ivp for a list of acceptable methods.')
            
        # Check CFL condition for diffusive component
        r = (self.D / self.c_p) * (dt / res**2)
        if r >= 1/2:
            warnings.warn('Unstable timestep for diffusion term.')

        # Initialize T_grid
        ts = dt*np.arange(0,N)
        phis = np.arange(-90, 90+res, res)
        J = len(phis)
        
        data = np.zeros((N, len(phis)))
        self.T_grid = xr.DataArray(data, dims=('time', 'latitude'), coords={'time': ts, 'latitude': phis})
        self.T_grid.loc[{'time': 0}] = [T_0(phi) for phi in phis]
        
        # Timestep T_grid
        for n in range(N-1):
            # Define previous T values
            T_last = self.T_grid.isel(time=n)

            # Define RHS of ODE (without diffusion)
            F = lambda t, T : 1/self.c_p * (self.S_0(phis) * (1 - self.alpha(T)) 
                                       - (1 - self.epsilon/2) * sigma * T**4)
            
            ode_sol = solve_ivp(F, (ts[n], ts[n+1]), np.array(T_last), method=method, t_eval=(ts[n], ts[n+1]))
            T_ode = ode_sol.y[:,1]

            # Diffusive component
            T_diffuse = xr.zeros_like(T_last)
            T_diffuse.loc[{'latitude':phis[0]}] = r*T_last.isel(latitude=1) + (1 - r)*T_last.isel(latitude=0)
            T_diffuse.loc[{'latitude':phis[-1]}] = (1 - r)*T_last.isel(latitude=-1) + r*T_last.isel(latitude=-2)
            for j in range(1, J-1):
                T_diffuse.loc[{'latitude':phis[j]}] = r*T_last.isel(latitude=j+1) + (1 - 2*r)*T_last.isel(latitude=j) + r*T_last.isel(latitude=j-1)
            T_diffuse = T_diffuse - T_last # Want to know the change from the present state

            # Superimpose diffusion plus ODE components to get next timestep:
            self.T_grid.loc[{'time': ts[n+1]}] = T_ode + T_diffuse
            
        return self.T_grid
        
    def show_animation(self):
        """
        Shows an animation of the simulation.
        """
        phis = self.T_grid['latitude']
        ts = self.T_grid['time']

        fig, ax = plt.subplots()
        ax.set(xlim=(-90, 90), ylim=(250, 300), 
            xlabel='Latitude', ylabel='Temperature')

        line = ax.plot(phis, self.T_grid.isel(time=0), color='g')[0]

        def animate(i):
            line.set_ydata(self.T_grid.isel(time=i))

        anim = FuncAnimation(fig, animate, interval=100, frames=len(ts)-1)
        
        plt.draw()
        plt.show()
        
        
        
        
        