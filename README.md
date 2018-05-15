# Vortex_Id - A 2D vortex Identification Code
===========================================

## Author: José Hugo Elsas
## Collaborator: Luca Moriconi 

Related publications: 
   - Elsas, J. H. and Moriconi. L. . Vortex identification from local properties of the vorticity field. Physics of Fluids 29, 015101 (2017);
   
---------------------------------------------------------------------------
# Content 

The objective of this code is to provide a way to systematically identify,
count and raise statistics of vortical structures in turbulent flow.

The input is a 2D velocity field, which can either be a true 2D flow or a 
cross cut from a 3D velocity field. As of the present moment, this flows include
synthetic 2D flows used to validate the identification part of the code and, 
ongoing, input from openFOAM files.

### Vortex_Id needs the following libraries and programs installed in your system

   - Cmake
   - GSL  : Gnu Scientific Library
   - inih : Ini file reader Library
   - numpy & matplotlib : python array and ploting libraries
   - A C99 compliant compiler (gcc 4.8.4 and icc 16.0.0 were tested)

### Instalation instructions:

   - Create folders named obj (for .o files), bin (for executables) and data (for results)
   - run make for main essay executable or main <executable> for the executable of your choice
   - If everything went well, everything should compile, even with some warnings it is not problematic.
   - All executables are stored on the bin folder

Obs:

    In case makefile don't work, try creating a obj folder on the home directory. It's necessary to store the .o files temporarely
