INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= seriesvtk
UTIL		= wholevtk
UTIL		= fftwtest
UTIL		= randIC

TRANSFORM	= fftw5
MODES		= M4
MODSOBJ		= mpi.o parameters.o \
		modes.o variables.o transform.o velocity.o \
		turb.o io.o 
#TURING COMPILATION
#COMPILER        = mpixlf90_r
#COMPFLAGS	= -qsuffix=cpp=f90 -c -O3 -ffree-line-length-none #-c -O3 -fbounds-check 
#LIBS		= -lm 
#ADA COMPILATION
COMPILER        = mpif90
COMPFLAGS       = -ffree-line-length-none -x f95-cpp-input -c -O3 -mcmodel=large \
	-I/network/aopp/cirrus/pred/sw/netcdf/4.3.3.1/gnu/492/include \
	-I/network/home/aopp/chantry/code/fftw_cirrus/include
LIBS            = -lm \
	-L/network/aopp/cirrus/pred/sw/netcdf/4.3.3.1/gnu/492/lib -lnetcdff -lnetcdf \
	-L/network/home/aopp/chantry/code/fftw_cirrus/lib -lfftw3f
#COMPILER        = mpiifort
#COMPFLAGS       = -fpp -132 -c -O3 -mcmodel=medium -I/network/software/ubuntu_trusty/netcdf/3.6.3/include
#LIBS            = -lm -L/network/software/ubuntu_trusty/netcdf/3.6.3/lib -lnetcdf
#------------------------------------------------------------------------

all : 	$(MODSOBJ) $(PROGDIR)main.f90
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)main.f90
	$(COMPILER) -o ./main.out main.o $(MODSOBJ) $(LIBS)

install : main.out
	if test ! -d $(INSTDIR); then mkdir -p $(INSTDIR); fi
	mv ./main.out $(INSTDIR)
	date > $(INSTDIR)/main.info
	echo $(HOSTNAME) >> $(INSTDIR)/main.info
	pwd >> $(INSTDIR)/main.info
	echo $(COMPILER) $(COMPFLAGS) >> $(INSTDIR)/main.info
	echo $(MODES) >> $(INSTDIR)/main.info
	grep define parallel.h | grep _Np >> $(INSTDIR)/main.info
	cut -d! -f1 $(PROGDIR)parameters.f90 | grep = | \
	   cut -d: -f3  >> $(INSTDIR)/main.info

util : 	$(MODSOBJ) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) $(COMPFLAGS) $(UTILDIR)/$(UTIL).f90
	$(COMPILER) -o ./$(UTIL).out $(UTIL).o $(MODSOBJ) $(LIBS)

#------------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.il core *.out

#------------------------------------------------------------------------
io.o : $(PROGDIR)io.f90 
	$(COMPILER) $(COMPFLAGS) -o io.o $(PROGDIR)io.f90

turb.o : $(PROGDIR)turb.f90 velocity.o mpi.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)turb.f90

mpi.o : $(PROGDIR)mpi.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)mpi.f90

parameters.o :	$(PROGDIR)parameters.f90 parallel.h
		$(COMPILER) $(COMPFLAGS) $(PROGDIR)parameters.f90

timestep.o : $(PROGDIR)timestep.f90 variables.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)timestep.f90

transform.o : $(PROGDIR)transform.$(TRANSFORM).f90 variables.o
	$(COMPILER) $(COMPFLAGS) -o transform.o \
	$(PROGDIR)transform.$(TRANSFORM).f90

variables.o : $(PROGDIR)variables.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) -o variables.o $(PROGDIR)variables.f90

velocity.o : $(PROGDIR)velocity.f90 transform.o modes.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)velocity.f90

modes.o : $(PROGDIR)modes.$(MODES).f90 parameters.o
	$(COMPILER) $(COMPFLAGS) -o modes.o \
	$(PROGDIR)modes.$(MODES).f90
