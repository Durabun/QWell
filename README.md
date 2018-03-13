# QWell
A quantum well sim

Being developed by Jonathan Eugenio and Mark Romero

If you have any questions, feel free to email me at: eugenejon1@gmail.com


Have you ever wanted to see the wavefunction for a particle in a one dimensional potential well?

Well now you can, with this app, coming soon!

## Current Features
* User defines the scenario, with the mass of the particle, the energy depth (in eV), and the half width well
* Plots each energy state with the representative wave function or density distribution function
* With Kivy, there is now a graphical interface for inputting parameters (9/6/2017)
* The Linear Potential Well has been completed (3/13/2018)
* The Multi-Barrier System has been completed (3/13/2018)

## In Development
* Linear Potential Well, with root finder and constants defined (11/29/2017)
* LPW wavefunction plots are currently in development
* Ideally this would be ported into the Kivy interface
* MBS is being developed for sloped barriers
* Looking into moving away from app development and into a flask web app form

## Future features
* Incorporate user defined perturbations
* More particles = more fun!

## Running the Simulation
Before starting, make sure you have the repo cloned onto your machine. There are also certain dependencies needed . These include matplotlib, numpy, scipy, kivy, and flask.

Make sure you are in QWell/

To run the Finite Square Well simulator in the terminal, type:
> $python quantumWell.py

To run the Finite Linear Well simulator in the terminal, type:
> $cd Airy/
> $python linPot.py

To run the Multi-Barrier System simulator in the terminal, type:
> $cd MBS/
> $python MBSgamma.py
Note that the MBS simulator is still being developed, hence the file name.

## Known Issues
There are a few issues with each part of the suite.
* The Finite Square Well simulator has a bug where the first and sometimes second wavefunctions are plotted for larger scenarios
* The Finite Linear Well simulator has a bug where the outer regions of the wavefunctions become discontinuous in scenarios with more than ~20 bound states.
* The Multi-Barrier System simulator cannot have any of the barrier heights equal to the particle energy, as this will cause a divide by zero error.
