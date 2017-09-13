# Kelp light model Makefile

##################
# Makefile links #
##################

# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# http://mrbook.org/blog/tutorials/make/
# http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules

###################
# Some Guidelines #
###################

# Every non-executing fortran file is a module
# use only the highest level modules needed
# Every module has an object
# Object rule dependencies should be identical to use statements in source
# Executable rule dependencies should be identical to use statements in source
# Executable compile line should include full tree of dependencies


#########
# Flags #
#########

# Project directories
BASE=.
BIN=$(BASE)/bin
SRC=$(BASE)/fortran
INC=$(BASE)/include
# External modules
EXT=$(SRC)/download

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

all: test old

pykelp3d: $(INC)/kelp3d.o $(INC)/hdf5_utils.o
	 $(H5FC) $(BFLAGS) $(SRC)/pykelp3d.f90 $(INC)/prob.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/kelp3d.o $(INC)/kelp_context.o $(INC)/hdf5_utils.o -o $(BIN)/pykelp3d


#########
# Tests #
#########

test: test_context test_gl test_grid test_gmres test_prob test_kelp_3d

$(INC)/test_context: $(INC)/utils.o $(INC)/kelp_context.o $(INC)/prob.o
	 $(FC) $(BFLAGS) $(SRC)/test_context.f90 $(INC)/prob.o $(INC)/utils.o $(INC)/kelp_context.o -o $(BIN)/test_context

test_gl: $(INC)/fastgl.o
	 $(FC) $(BFLAGS) $(SRC)/test_gl.f90 $(INC)/fastgl.o -o $(BIN)/test_gl

test_grid: $(INC)/kelp_context.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/prob.o
	 $(FC) $(BFLAGS) $(SRC)/test_grid.f90 $(INC)/kelp_context.o $(INC)/fastgl.o $(INC)/prob.o $(INC)/sag.o $(INC)/utils.o -o $(BIN)/test_grid

test_gmres: $(INC)/mgmres.o $(INC)/hdf5_utils.o
	 $(H5FC) $(BFLAGS) $(SRC)/test_gmres.f90 $(INC)/mgmres.o $(INC)/utils.o $(INC)/hdf5_utils.o $(INC)/kelp_context.o $(INC)/sag.o $(INC)/fastgl.o $(INC)/prob.o -o $(BIN)/test_gmres 

$(INC)/test_prob: $(INC)/prob.o
	 $(FC) $(BFLAGS) $(SRC)/test_prob.f90 $(INC)/prob.o -o $(BIN)/test_prob

test_kelp3d: $(INC)/test_kelp3d_mod.o
	 $(H5FC) $(BFLAGS) $(SRC)/test_kelp3d.f90 $(INC)/test_kelp3d_mod.o $(INC)/prob.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/kelp3d.o $(INC)/kelp_context.o -o $(BIN)/test_kelp3d

test_rte3d: $(INC)/test_kelp3d_mod.o $(INC)/rte3d.o
	 $(H5FC) $(BFLAGS) $(SRC)/test_rte3d.f90 $(INC)/rte_sparse_matrices.o $(INC)/mgmres.o $(INC)/rte3d.o $(INC)/test_kelp3d_mod.o $(INC)/prob.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/kelp3d.o $(INC)/kelp_context.o -o $(BIN)/test_rte3d

### Old tests
old: test_interp test_rte2d test_vsf

test_interp: $(INC)/utils.o
	 $(FC) $(BFLAGS) $(SRC)/test_interp.f90 $(INC)/utils.o -o $(BIN)/test_interp

test_rte2d: $(INC)/rte2d.o $(INC)/rte_core.o $(INC)/utils.o
	 $(FC) $(BFLAGS) $(SRC)/test_rte2d.f90 $(INC)/rte2d.o $(INC)/rte_core.o $(INC)/utils.o -o $(BIN)/test_rte2d

test_vsf: $(INC)/rte_core.o $(INC)/utils.o
	 $(FC) $(BFLAGS) $(SRC)/test_vsf.f90 $(INC)/rte_core.o $(INC)/utils.o -o $(BIN)/test_vsf


################
# Object files #
################

$(INC)/sag.o: $(INC)/utils.o $(INC)/fastgl.o
	$(FC) $(OFLAGS) $(SRC)/sag.f90 -o $(INC)/sag.o

$(INC)/hdf5_utils.o: $(INC)/utils.o $(INC)/kelp_context.o
	$(H5FC) $(OFLAGS) $(SRC)/hdf5_utils.f90 -o $(INC)/hdf5_utils.o

$(INC)/kelp_context.o: $(INC)/sag.o $(INC)/prob.o
	$(FC) $(OFLAGS) $(SRC)/kelp_context.f90 -o $(INC)/kelp_context.o

$(INC)/kelp3d.o: $(INC)/kelp_context.o
	$(FC) $(OFLAGS) $(SRC)/kelp3d.f90 -o $(INC)/kelp3d.o

$(INC)/test_kelp3d_mod.o: $(INC)/kelp3d.o $(INC)/hdf5_utils.o
	$(FC) $(OFLAGS) $(SRC)/test_kelp3d_mod.f90 -o $(INC)/test_kelp3d_mod.o

$(INC)/rte_sparse_matrices.o: $(INC)/sag.o $(INC)/kelp_context.o $(INC)/mgmres.o
	$(FC) $(OFLAGS) $(SRC)/rte_sparse_matrices.f90 -o $(INC)/rte_sparse_matrices.o

$(INC)/rte3d.o: $(INC)/kelp_context.o $(INC)/rte_sparse_matrices.o
	$(FC) $(OFLAGS) $(SRC)/rte3d.f90 -o $(INC)/rte3d.o

# Old

$(INC)/rte_core.o: $(INC)/utils.o
	$(FC) $(OFLAGS) $(SRC)/rte_core.f90 -o $(INC)/rte_core.o

$(INC)/utils.o:
	$(FC) $(OFLAGS) $(SRC)/utils.f90 -o $(INC)/utils.o

$(INC)/rte2d.o: $(INC)/rte_core.o
	$(FC) $(OFLAGS) $(SRC)/rte2d.f90 -o $(INC)/rte2d.o

# External

$(INC)/fastgl.o:
	$(FC) $(OFLAGS) $(EXT)/fastgl.f90 -o $(INC)/fastgl.o

$(INC)/mgmres.o:
	$(FC) $(OFLAGS) $(EXT)/mgmres.f90 -o $(INC)/mgmres.o

$(INC)/prob.o:
	$(FC) $(OFLAGS) $(EXT)/prob.f90 -o $(INC)/prob.o

# Utils

clean: rmo
	rm -f $(INC)/*.mod $(INC)/*.o $(BIN)/*

rmo: 
	rm -f $(BASE)/*.o $(BASE)/*.mod

ls:
	ls $(SRC) $(BIN) $(INC)
