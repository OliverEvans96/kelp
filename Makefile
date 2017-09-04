# Kelp light model Makefile

##################
# Makefile links #
##################

# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# http://mrbook.org/blog/tutorials/make/
# http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules

#########
# Flags #
#########

# Project directories
BASE=.
BIN=$(BASE)/bin
SRC=$(BASE)/fortran
INC=$(BASE)/include

# Fortran Compiler
FC=gfortran
# HDF5 Fortan Compiler
H5FC=h5fc

# Fortran Compilation flags
# Object files (.o)
OFLAGS=-J$(INC) -I$(INC) -c
# Binary files (executable)
BFLAGS=-J$(INC) -I$(INC)

###############
# Executables #
###############

all: test_interp test_rte2d test_vsf test_context

test_interp: utils.o
	 $(FC) $(BFLAGS) $(SRC)/test_interp.f90 $(INC)/utils.o -o $(BIN)/test_interp

test_context: utils.o kelp_context.o
	 $(FC) $(BFLAGS) $(SRC)/test_context.f90 $(INC)/utils.o $(INC)/kelp_context.o -o $(BIN)/test_context

test_gl: 
	 $(FC) $(BFLAGS) $(SRC)/test_gl.f90 $(SRC)/download/fastgl.f90 -o $(BIN)/test_gl

test_gmres: hdf5_utils.o utils.o
	 $(H5FC) $(BFLAGS) $(SRC)/test_gmres.f90 $(SRC)/download/mgmres.f90 $(INC)/utils.o $(INC)/hdf5_utils.o -o $(BIN)/test_gmres 

test_prob: 
	 $(FC) $(BFLAGS) $(SRC)/test_prob.f90 $(SRC)/download/prob.f90 -o $(BIN)/test_prob

test_rte2d: rte2d.o
	 $(FC) $(BFLAGS) $(SRC)/test_rte2d.f90 $(INC)/rte2d.o $(INC)/rte_core.o $(INC)/utils.o -o $(BIN)/test_rte2d

test_vsf: rte_core.o
	 $(FC) $(BFLAGS) $(SRC)/test_vsf.f90 $(INC)/rte_core.o $(INC)/utils.o -o $(BIN)/test_vsf

################
# Object files #
################

rte2d.o: rte_core.o
	$(FC) $(OFLAGS) $(SRC)/rte2d.f90 -o $(INC)/rte2d.o

sag.o:
	$(FC) $(OFLAGS) $(SRC)/sag.f90 -o $(INC)/sag.o

hdf5_utils.o: utils.o
	$(H5FC) $(OFLAGS) $(SRC)/hdf5_utils.f90 -o $(INC)/hdf5_utils.o

kelp_context.o: sag.o
	$(FC) $(OFLAGS) $(SRC)/kelp_context.f90 -o $(INC)/kelp_context.o

rte_core.o: utils.o
	$(FC) $(OFLAGS) $(SRC)/rte_core.f90 -o $(INC)/rte_core.o

utils.o:
	$(FC) $(OFLAGS) $(SRC)/utils.f90 -o $(INC)/utils.o

clean: rmo
	rm -f $(INC)/*.mod $(INC)/*.o $(BIN)/*

rmo: 
	rm -f $(BASE)/*.o


ls:
	ls $(SRC) $(BIN) $(INC)
