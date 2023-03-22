# Computational Fluid Dynamics - Advection
This program solves an advection problem using runge-kutta methods for integration of the CDF ODEs, for a gaussian or top hat profile using diferent slope delimiters. 

## Instructions
Advection module:
* **.Grid1()** generates the 1D grid and initial set up: ng (gosths cells), nx (number of cells) and x limits.
* **.Simulation()** takes an object as input
* **initial_conditions()** will fill the first array of or grid with initial conditions using a "top hat" or "gaussian" profile.
* **time_loop()** will perform the operations iteratively up to a max-time.

#### Example: Top hat profile using the slope type "minmod"

- Initialize **Grid1() if no x limits are given, they are set up to 0 and 1. A grid is created within the class.

- Then make the grid boundaries periodic by calling the function fill_bcs()

- We call .Simulation() providing the object return from .Grid1(), the cfl constant c, u, and the slope type which in this case we use "minmod".

- We apply the initial conditions on the revious object by using initial_conditions() which acepts, "tophat" or "gaussian"

- Finally we loop over a t_max establish by defaut as the period * number of periods. This function will operate over the grid "a" that will store the results.

- Plot g.x (x axis) vs g.a (results)

## Instalations
- Download the folder "cfdadvection"
- Install the dependencies
- Import cdfadvection.py

## Contact
For questions or feedback, please contact me at [brandon.minta@yachaytech.edu.ec].
