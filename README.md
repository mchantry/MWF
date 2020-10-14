This code simulated Model Waleffe flow in a retangular domain.

To build:
1. Edit Makefile to point to netcdf & fftw
2. make;make install

To generate random initial condition:
1. Ensure UTIL in makefile is set as randIC
2. make util
3. ./randIC.out
4. This produces state0000.cdf.dat
5. Rename to state.cdf.in for use as initial condition.
6. If state relaminarises then increase scl & recompile, or run at higher Re.

To run.
1. Copy main.info, main.out & state.cdf.in to folder.
2. Run ./main.out
3. vel_energy.dat keeps a running output of the total energy.
3. Delete RUNNING file to softly kill run.
4. For parallel simulation edit _Np in parallel.h & recompile.

To control.
1. program/parameters.f90 contains all parameters to be edited.
2. Recompile.
3. main.info has a copied of the parameters used at compilation.

To plot output.
1. In matlab [x,z,u]=GridUy('state.cdf.in','U',0.) extracts the U field at y=0.
2. contourf(x,z,u) to visualise.

Questions?
Feel free to contact me.