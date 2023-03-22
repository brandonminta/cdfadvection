# Computational Fluid Dynamics - Advection
This program solves an advection problem using runge-kutta methods for integration of the CDF ODEs, for a gaussian or top hat profile using diferent slope delimiters. 

## Instructions
Advection module:
* **.Grid1()** generates the 1D grid and initial set up: ng (gosths cells), nx (number of cells) and x limits.
* **.Simulation()** takes an object as input
* **initial_conditions()** will fill the first array of or grid with initial conditions using a "top hat" or "gaussian" profile.
* **time_loop()** will perform the operations iteratively up to a max-time.



## Instalations
- Download the folder "cfdadvection"
- Install the dependencies: pip install -r requirements.txt
- Import cdfadvection.py

## Contact
For questions or feedback, please contact me at [brandon.minta@yachaytech.edu.ec].
