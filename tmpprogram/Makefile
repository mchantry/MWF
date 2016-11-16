INSTDIR		= ./install/
PROGDIR		= ./program/
UTILDIR		= ./utils/
UTIL		= seriesvtk
#UTIL		= wholevtk
#UTIL		= turbfrac
#UTIL		= interpstate

TRANSFORM	= fftw2
MODES		= M4

MODSOBJ		= io.o mpi.o modes.o parameters.o \
		  	variables.o transform.o velocity.o 
#MODSOBJ         = io.o meshs.o mpi.o parameters.o \
#                  timestep.o transform.o variables.o velocity.o

#COMPILER        = /opt/openmpi/bin/mpif77
#COMPILER	= mpif77
COMPILER	= mpiifort
#COMPILER        = /opt/gcc-4.6.3/bin/gfortran
COMPFLAGS	= -132 -cpp -c -O3 -mcmodel=medium
#COMPFLAGS	= -ffree-line-length-none -static -x f95-cpp-input -c -O3 \
			 -I/usr/local/Cellar/netcdf/4.3.3.1_4/include 
                  #-C #-pg
LIBS		= -lm
#LIBS		= -lfftw3 -lnetcdff -lnetcdf \
		-L/usr/local/Cellar/fftw/3.3.4_1/lib/ \
		-L/usr/local/Cellar/netcdf/4.3.3.1_4/lib/

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
io.o : $(PROGDIR)io.f90 velocity.o 
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)io.f90

meshs.o : $(PROGDIR)meshs.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)cheby.f
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)meshs.f90

mpi.o : $(PROGDIR)mpi.f90 parallel.h
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)mpi.f90

parameters.o :	$(PROGDIR)parameters.f90 parallel.h
		$(COMPILER) $(COMPFLAGS) $(PROGDIR)parameters.f90

timestep.o : $(PROGDIR)timestep.f90 variables.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)timestep.f90

transform.o : $(PROGDIR)transform.$(TRANSFORM).f90 variables.o
	$(COMPILER) $(COMPFLAGS) -o transform.o \
	$(PROGDIR)transform.$(TRANSFORM).f90

modes.o : $(PROGDIR)modes.$(MODES).f90 parameters.o
	$(COMPILER) $(COMPFLAGS) -o modes.o \
	$(PROGDIR)modes.$(MODES).f90

variables.o : $(PROGDIR)variables.f90 parameters.o mpi.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)variables.f90

velocity.o : $(PROGDIR)velocity.f90 transform.o modes.o
	$(COMPILER) $(COMPFLAGS) $(PROGDIR)velocity.f90

