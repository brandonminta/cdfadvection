import numpy as np


def maxmod(a, b):
    if abs(a) > abs(b) and a*b > 0.:
        return a
    elif abs(b) > abs(a) and a*b > 0.:
        return b
    else:
        return 0.
    
def minmod(a, b):
    if abs(a) < abs(b) and a*b > 0.:
        return a
    elif abs(b) < abs(a) and a*b > 0.:
        return b
    else:
        return 0.    


class Grid1(object):
    
    
    def __init__(self, ng, nx, xmin = 0., xmax = 1.):
        self.ng = ng               #number of ghost cells on each side
        self.nx = nx               #number of cells
        self.xmin = xmin           #limits
        self.xmax = xmax           #limits
        
        # Lowest and highest indices
        self.ilo = ng
        self.ihi = ng + nx - 1
        
        #step size
        self.dx = (xmax - xmin)/(nx)
        
        # Empty vector for the solution in double precision:
        self.a = np.zeros((nx + 2*ng), dtype = np.float64)
        
        # X axis vector -> we need the cell centres:
        self.x = xmin + (np.arange(nx + 2*ng) - ng + 0.5)*self.dx
        
        
    #function to Fill the boundary conditions
    def fill_bcs(self):
        # Fill up the boundary conditions BCs
        for n in range(self.ng):
            # The left BC is
            self.a[self.ilo - 1 - n] = self.a[self.ihi - n]
            # The right BC is
            self.a[self.ihi + 1 + n] = self.a[self.ilo + n]
            
    def norm(self, value):
        if len(value) != (2*self.ng + self.nx):
            return None

        return np.sqrt(self.dx*np.sum(value[self.ilo:self.ihi+1]**2)) # L2 norm

    
    
class Simulation(object):
    
    def __init__(self, grid, u, c = 0.7, slope_type = "centered"):
        
        self.grid = grid
        self.t = 0.
        self.u = u           # Define perturbation velocity
        self.c = c           # Define Courant number
        self.slope_type = slope_type
    
    # Get time step
    def timestep(self):
        return self.c*self.grid.dx/self.u
    
    # Calculate period
    def period(self):
        return (self.grid.xmax - self.grid.xmin)/self.u
    
    
    def initial_conditions(self, type="tophat"):
        if type == "tophat":
            self.grid.a[:] = 0.0
            self.grid.a[np.logical_and(self.grid.x >= 1./3., self.grid.x <= 2./3.)] = 1.0

        elif type == "gaussian":
            #gaussian with mean 0.5 and standar deviation 0.1
            self.grid.a[:] = np.exp(-(self.grid.x - 0.5)**2 / (2 * (0.1)**2)) / (0.1 * np.sqrt(2 * np.pi))

    
    def states(self, dt):

        # Empty vector for the slope
        slope = np.zeros((self.grid.nx + 2*self.grid.ng), dtype = np.float64)
        
            
        if self.slope_type == "godunov":
            # 1st approach: Godunov approach
            slope[:] = 0.0

        elif self.slope_type == "unlimited":
            # 2nd approach: Centred difference method
            for i in range(self.grid.ilo-1, self.grid.ihi+2):
                slope[i] = 0.5*(self.grid.a[i+1] - self.grid.a[i-1])/self.grid.dx

        elif self.slope_type == "minmod":
            # 3rd approach: minmod limiter
            for i in range(self.grid.ilo-1, self.grid.ihi+2):
                slope[i] = minmod( (self.grid.a[i] - self.grid.a[i-1])/self.grid.dx, 
                                  (self.grid.a[i+1] - self.grid.a[i])/self.grid.dx )

        elif self.slope_type == "MC":
            # 4th approach: MC limiter
            for i in range(self.grid.ilo-1, self.grid.ihi+2):
                slope[i] = minmod(minmod( 2.0*(self.grid.a[i] - self.grid.a[i-1])/self.grid.dx, 
                                         2.0*(self.grid.a[i+1] - self.grid.a[i])/self.grid.dx ),
                                  0.5*(self.grid.a[i+1] - self.grid.a[i-1])/self.grid.dx)

        elif self.slope_type == "superbee":
            # 5th approach: MC limiter
            for i in range(self.grid.ilo-1, self.grid.ihi+2):
                min1 = minmod( (self.grid.a[i+1] - self.grid.a[i])/self.grid.dx, 
                              2.0*(self.grid.a[i] - self.grid.a[i-1])/self.grid.dx )

                min2 = minmod( (self.grid.a[i] - self.grid.a[i-1])/self.grid.dx, 
                              2.0*(self.grid.a[i+1] - self.grid.a[i])/self.grid.dx )

                slope[i] = maxmod(min1, min2)
        
        elif self.slope_type == "vanleer":
            # 6th approach: Van Leer
            for i in range(self.grid.ilo-1, self.grid.ihi+2):
                solpe[i] = (((self.grid.a[i] - self.grid.a[i-1])/self.grid.dx) + \
                    abs((self.grid.a[i] - self.grid.a[i-1])/self.grid.dx))/(1+abs((self.grid.a[i] - self.grid.a[i-1])/self.grid.dx))
                
            
        # Empty vector for the L and R states
        al = np.zeros((self.grid.nx + 2*self.grid.ng), dtype = np.float64)
        ar = np.zeros((self.grid.nx + 2*self.grid.ng), dtype = np.float64)
        
        # Add derivative calculation
        for i in range(self.grid.ilo, self.grid.ihi+2):

            # Compute L state
            al[i] = self.grid.a[i-1] + 0.5*self.grid.dx*(1. - self.u*dt/self.grid.dx)*slope[i-1]

            # Compute R state
            ar[i] = self.grid.a[i] - 0.5*self.grid.dx*(1. + self.u*dt/self.grid.dx)*slope[i]

        return al, ar
        
    def riemann(self, al, ar):
        # Computing the flux
        if self.u > 0.:
            return self.u*al
        else:
            return self.u*ar
 

    def conservative_update(self, dt, flux):

        # Empty vector for updated solution
        anew = np.zeros((self.grid.nx + 2*self.grid.ng), dtype = np.float64)
        
        # Update
        anew[self.grid.ilo:self.grid.ihi+1] = self.grid.a[self.grid.ilo:self.grid.ihi+1] - \
            dt/self.grid.dx * (flux[self.grid.ilo +1 :self.grid.ihi + 2] - flux[self.grid.ilo:self.grid.ihi+1])

        return anew


    def time_loop(self, n_periods=5):
        self.t = 0.

        tmax = n_periods*self.period()


        # main evolution loop
        while self.t < tmax:

            # fill the boundary conditions
            self.grid.fill_bcs()

            # get the timestep
            dt = self.timestep()

            if self.t + dt > tmax:
                dt = tmax - self.t

            # get the interface states
            al, ar = self.states(dt)

            # solve the Riemann problem at all interfaces
            flux = self.riemann(al, ar)

            # do the conservative update
            anew = self.conservative_update(dt, flux)
            
            # Retrieve solution
            self.grid.a[:] = anew[:]
            
            # Advance in time
            self.t += dt

        
        
        