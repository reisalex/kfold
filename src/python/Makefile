# Makefile for KFOLD FORTRAN 2 PYTHON EXTENSION MODULE
#
# Authors: Alex Reis
# Copyright (c) 2016 (Please Refer to LICENCE)

F90 = gfortran
FFLAGS = -fPIC -O3 -fcheck=all

# ---------------------------------------------------------------------
# Debugging options
# ---------------------------------------------------------------------

# -g :: Produce debugging information in the operating system's native 
# format (stabs, COFF, XCOFF, or DWARF). GDB can work with this 
# debugging information.

# -Wall :: Enables commonly used warning options pertaining to usage
# that we recommend avoiding and that we believe are easy to avoid.

# -fcheck=all :: Enable the generation of run-time checks; the argument
# shall be a comma-delimited list of the following keywords.

# ---------------------------------------------------------------------
# Optimize options
# ---------------------------------------------------------------------

# -O3 :: Turns on all optimizations specified by -O2 and also turns on 
# the -finline-functions, -funswitch-loops, -fpredictive-commoning, 
# -fgcse-after-reload, -ftree-vectorize and -fipa-cp-clone options.

# ---------------------------------------------------------------------
# Profiling options
# ---------------------------------------------------------------------

# -pg :: Generate extra code to write profile information suitable for
# the analysis program gprof. You must use this option when compiling the 
# source files you want data about, and you must also use it when linking.

# ---------------------------------------------------------------------
# OBJECTS
# ---------------------------------------------------------------------

MAIN = run.o kfold.o

MODULES = rnavar.o

CLASS = class_rnafold.o

SUBS = ssareaction.o convert.o v2ct.o setupnuc.o readdata.o random.o

ENERGY = deltag_he.o deltag_hr.o deltag_hi.o deltag_hm.o deltag_hd.o \
         estruc.o eloop.o ebulge.o ehair.o estack.o edangle.o \
         tint11.o tint12.o tint22.o tloop.o tstack.o tstackh.o \
         tstacki.o tdangle5.o tdangle3.o

OBJECTS = $(MAIN) $(MODULES) $(CLASS) $(CLASS_KFOLD) $(SUBS) $(ENERGY)

# ---------------------------------------------------------------------
# MAKE COMMANDS
# ---------------------------------------------------------------------

all: clean $(OBJECTS)

clean:
	rm -f *.o *.mod

veryclean:
	rm -f *.o *.mod *.x

# ---------------------------------------------------------------------
# COMPILE COMMANDS
# ---------------------------------------------------------------------

# --------- MAIN ------------

run.o: run.f90 $(MODULES)
	$(F90) $(FFLAGS) -c run.f90

kfold.o: kfold.f90 $(MODULES) $(CLASS)
	$(F90) $(FFLAGS) -c kfold.f90

# --------- MODULES -----------

rnavar.o: ../rnavar.f90
	$(F90) $(FFLAGS) -c ../rnavar.f90

# --------- CLASSES -----------

class_rnafold.o: $(MODULES) ../class_rnafold.f90
	$(F90) $(FFLAGS) -c ../class_rnafold.f90

# --------- SUBS ------------
ssareaction.o: $(MODULES) $(CLASS) ssareaction.f90
	$(F90) $(FFLAGS) -c ssareaction.f90
convert.o: ../convert.f90
	$(F90) $(FFLAGS) -c ../convert.f90
v2ct.o: ../v2ct.f90
	$(F90) $(FFLAGS) -c ../v2ct.f90
setupnuc.o: $(MODULES) ../setupnuc.f90
	$(F90) $(FFLAGS) -c ../setupnuc.f90
readdata.o: $(MODULES) ../readdata.f90
	$(F90) $(FFLAGS) -c ../readdata.f90
random.o: random.f90
	$(F90) $(FFLAGS) -c random.f90

# --------- ENERGY ----------

deltag_he.o: $(MODULES) ../ENERGY/deltag_he.f90
	$(F90) $(FFLAGS) -c ../ENERGY/deltag_he.f90
deltag_hr.o: $(MODULES) ../ENERGY/deltag_hr.f90
	$(F90) $(FFLAGS) -c ../ENERGY/deltag_hr.f90
deltag_hi.o: $(MODULES) ../ENERGY/deltag_hi.f90
	$(F90) $(FFLAGS) -c ../ENERGY/deltag_hi.f90
deltag_hm.o: $(MODULES) ../ENERGY/deltag_hm.f90
	$(F90) $(FFLAGS) -c ../ENERGY/deltag_hm.f90
deltag_hd.o: $(MODULES) ../ENERGY/deltag_hd.f90
	$(F90) $(FFLAGS) -c ../ENERGY/deltag_hd.f90
estruc.o: $(MODULES) ../ENERGY/estruc.f90
	$(F90) $(FFLAGS) -c ../ENERGY/estruc.f90
eloop.o: $(MODULES) ../ENERGY/eloop.f90
	$(F90) $(FFLAGS) -c ../ENERGY/eloop.f90
ebulge.o: $(MODULES) ../ENERGY/ebulge.f90
	$(F90) $(FFLAGS) -c ../ENERGY/ebulge.f90
ehair.o: $(MODULES) ../ENERGY/ehair.f90
	$(F90) $(FFLAGS) -c ../ENERGY/ehair.f90
estack.o: $(MODULES) ../ENERGY/estack.f90
	$(F90) $(FFLAGS) -c ../ENERGY/estack.f90
edangle.o: $(MODULES) ../ENERGY/edangle.f90
	$(F90) $(FFLAGS) -c ../ENERGY/edangle.f90
tint11.o: ../ENERGY/tint11.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tint11.f90
tint12.o: ../ENERGY/tint12.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tint12.f90
tint22.o: ../ENERGY/tint22.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tint22.f90
tloop.o: ../ENERGY/tloop.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tloop.f90
tstack.o: ../ENERGY/tstack.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tstack.f90
tstackh.o: ../ENERGY/tstackh.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tstackh.f90
tstacki.o: ../ENERGY/tstacki.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tstacki.f90
tdangle5.o: ../ENERGY/tdangle5.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tdangle5.f90
tdangle3.o: ../ENERGY/tdangle3.f90
	$(F90) $(FFLAGS) -c ../ENERGY/tdangle3.f90
