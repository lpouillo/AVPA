PROG =	AV.exe

SRCS =	anisotropic_viscosity.f90 AV_calc.f90 AV_inout.f90 AV_variables.f90 decmod.f90 \
	nrmod.f90 tensmod.f90

OBJS =	anisotropic_viscosity.o AV_calc.o AV_inout.o AV_variables.o decmod.o nrmod.o \
	tensmod.o

LIBS =	

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = 
F90 = gfortran
F90FLAGS = 
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

anisotropic_viscosity.o: AV_calc.o AV_inout.o AV_variables.o tensmod.o
AV_calc.o: AV_inout.o AV_variables.o decmod.o tensmod.o
AV_inout.o: AV_variables.o 
decmod.o: nrmod.o tensmod.o
