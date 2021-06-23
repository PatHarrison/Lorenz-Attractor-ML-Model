"""
Creates a Lorenz attractor (butterfly) object
"""
import numpy as np
import scipy.integrate

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button

import warnings

class LorenzAttractor:
    """ Lorenz Attractor object """

    def __init__(self, t, initialState, sigma, beta, rho):
        """
        Constructor Method for LorenzAttractor object. Class of methods for
        the computation of

            dX/dt = \sigma(Y-X)
            dY/dt = -XZ+\rho X-Y
            dZ/dt = XY-\beta Z

        Args:
            initialState: 3-tuple, (X, Y, Z)=Initial state of the attractor
            sigma: float, represents physical parameter $\sigma$
            beta: float, represents physical parameter $\beta$
            rho: float, represents physical paraneter $\rho$

        Returns:
            None
        """

        self._initialState = initialState
        self._sigma = sigma
        self._beta = beta
        self._rho = rho
        
        self._t = t
        
        self._states = self.__solve(t)
        
        self.update_critical_points()

        return None
    
    def __repr__(self):
        return f""" Lorenz Attractor: (s={self._sigma}, b={self._beta}, 
                    r={self._rho}), L0={self._initialState}
                """

    def __derivatives(self, state, t0):
        """ Puts system of equations into a length 3 array (Xdot, Ydot, Zdot).
            t0 is needed for the integration step implimintaion
        """

        x, y, z = state
        return np.array([self._sigma*(y-x),
                         x*(self._rho-z)-y,
                         x*y-self._beta*z])

    def __solve(self, t):
        """ Integrates to get state of system in time steps specified by t

        Args:
            t: A sequence of time points to solve for the state.
        """

        return scipy.integrate.odeint(self.__derivatives,
                                      self._initialState,
                                      t)

    def set_params(self, sigma, beta, rho):
        """ Setter function for the physical paramters and updates 
            states attributes
        """
        
        self._sigma = sigma
        self._beta = beta
        self._rho = rho
        
        self._states = self.__solve(self._t)
    
    def get_initialState(self):
        return self._initialState
        
    def get_sigma(self):
        return self._sigma
    
    def get_beta(self):
        return self._beta
    
    def get_rho(self):
        return self._rho
        
    def get_states(self):
        """ Getter method for states 
        Returns:
            X, Y, Z
        """
        return self._states[:,0], self._states[:,1], self._states[:,2]
    
    def get_critical_points(self):
        """ Shape (2, 3): [(X,Y,Z), (X,Y,Z)]"""
        self.update_critical_points()
        return self._critical_points
    
    def update_critical_points(self):
        """ Calcualte critical points if self._rho > 1 """
        if self.check_stable_critical_points():
            #Calculate critical points (self._rho > 1)
            self._critical_points = np.array([(np.sqrt(self._beta*(self._rho-1)),
                                      np.sqrt(self._beta*(self._rho-1)),
                                      self._rho-1),
                                     (-np.sqrt(self._beta*(self._rho-1)),
                                      -np.sqrt(self._beta*(self._rho-1)),
                                      self._rho-1)
                                    ])
        else:
            self._critical_points = np.empty((2,3)) # None type objects

    def check_stable_critical_points(self):
        """ get critical points
            https://en.wikipedia.org/wiki/Lorenz_system
        """
        if self._rho > 1:
            return True
        else:
            warnings.warn("No stable critical points")
            return False
        
    def plot(self):
        """ Plots the states as components and in 3-space
        
        returns:
            plot object
        """
        
        fig = plt.figure(constrained_layout=True, figsize=plt.figaspect(.75))
        fig.tight_layout()
        gs = GridSpec(6, 2, 
                      figure=fig,
                      right=1,
                      height_ratios=[1,1,1,.5,.5,.5])

        title = r"Lorenz attractor: $\sigma=${}, $\beta=${}, $\rho=${}"
        
        # Set up subplot layout
        ax1 = fig.add_subplot(gs[0,0], ylabel="X")
        ax2 = fig.add_subplot(gs[1,0], ylabel="Y")
        ax3 = fig.add_subplot(gs[2,0], ylabel="Z")
        ax4 = fig.add_subplot(gs[:,1], projection='3d',
                             xlabel="X", ylabel="Y", zlabel="Z")
        
        axes = (ax1, ax2, ax3, ax4)
        axes[2].set_xlabel('Time')

        # Slider axes
        axSigma = fig.add_subplot(gs[3,0])
        axBeta = fig.add_subplot(gs[4,0])
        axRho = fig.add_subplot(gs[5,0])
        
        # create sliders
        sigmaSlider = Slider(axSigma, 
                             label=r"$\sigma$", 
                             valmin=0, 
                             valmax=20, 
                             valinit=10, 
                             valstep=0.5)
        betaSlider = Slider(axBeta, 
                             label=r"$\beta$", 
                             valmin=0, 
                             valmax=10, 
                             valinit=8/3.0, 
                             valstep=0.1)
        rhoSlider = Slider(axRho, 
                             label=r"$\rho$", 
                             valmin=0, 
                             valmax=50, 
                             valinit=28, 
                             valstep=0.5)
        
        def update(val):
            """ update function for sliders. Plots the new states
                CAUTION: will change the states attribute and critical points.

                removing old lines is funky and critical point markers did not 
                get removed hence the redundant removal of lines.
                
               returns:
                   figure, ax, and tuple of sliders. Otherwise matplotlib 
                   will not update with sliders inside a function.
                   https://github.com/matplotlib/matplotlib/issues/3105
            """
            self.set_params(sigmaSlider.val, betaSlider.val, rhoSlider.val)
            self.update_critical_points()
            cPoint1, cPoint2 = self.get_critical_points()

            for indx in range(4):
                [l.remove() for l in axes[indx].lines] # remove old lines
            
            for indx in range(3): # plot new Lines
                [l.remove() for l in axes[indx].lines] # make sure it's removed
                axes[indx].plot(self._t, self._states[:,indx], 
                                lw=0.75, color='k')

            # plot 3d projection of butterfly
            [l.remove() for l in axes[3].lines] # make sure removed
            axes[3].plot(self._states[:,0],
                         self._states[:,1],
                         self._states[:,2],
                         color='k', lw=0.75)
            axes[3].plot(cPoint1[0],cPoint1[1],cPoint1[2],"ro-")
            axes[3].plot(cPoint2[0],cPoint2[1],cPoint2[2],"ro-")

            # update title
            fig.suptitle(title.format(str(self._sigma), 
                                      str(self._beta), 
                                      str(self._rho)),
                        fontsize=16)
        update(1)
        sigmaSlider.on_changed(update)
        betaSlider.on_changed(update)
        rhoSlider.on_changed(update)
        
        #show plot
        plt.show()

        return fig, axes, (sigmaSlider, betaSlider, rhoSlider)


if __name__=="__main__":
    print("Starting Lorenz Object")
    tStart = 0.0
    tEnd = 30.0
    dt = 0.003

    t = np.arange(tStart, tEnd, dt)
    
    LA = LorenzAttractor(t, (1.0,1.0,1.0))
    print("plot")
    LA.plot()
