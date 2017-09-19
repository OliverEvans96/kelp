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

# Profiling (timing) flags
# Add to OFLAGS and BFLAGS
# (Should replace CFLAGS in OFLAGS)
# Compile & run normally,
# Then call gprof `executable`
PFLAGS=-g -O0

# Optimize performance
#https://stackoverflow.com/questions/42386065/inlining-functions-in-fortran
#OPTFLAGS=-Ofast -flto

# Normal compilation flag (no profiling)
CFLAGS=-c

# Fortran Compilation flags
# Object files (.o)
OFLAGS=-J$(INC) -I$(INC) $(CFLAGS) $(PFLAGS) $(OPTFLAGS)
# Binary files (executable)
BFLAGS=-J$(INC) -I$(INC) $(PFLAGS) $(OPTFLAGS)


###############
# Executables #
###############

all: test old

pykelp3d: $(SRC)/pykelp3d.f90 $(INC)/prob.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/kelp3d.o $(INC)/kelp_context.o $(INC)/hdf5_utils.o
	$(H5FC) $(BFLAGS) $^ -o $@

#########
# Tests #
#########

test: test_context test_gl test_grid test_gmres test_prob test_kelp_3d

test_context: $(SRC)/test_context.f90 $(INC)/prob.o $(INC)/utils.o $(INC)/kelp_context.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_gl: $(SRC)/test_gl.f90 $(INC)/fastgl.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_grid: $(SRC)/test_grid.f90 $(INC)/kelp_context.o $(INC)/fastgl.o $(INC)/prob.o $(INC)/sag.o $(INC)/utils.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_gmres: $(SRC)/test_gmres.f90 $(INC)/mgmres.o $(INC)/utils.o $(INC)/hdf5_utils.o $(INC)/kelp_context.o $(INC)/sag.o $(INC)/fastgl.o $(INC)/prob.o
	$(H5FC) $(BFLAGS) $^ -o $(BIN)/$@
$(INC)/test_prob: $(SRC)/test_prob.f90 $(INC)/prob.o $(INC)/prob.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_kelp3d: $(SRC)/test_kelp3d.f90 $(INC)/test_kelp3d_mod.o $(INC)/prob.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/kelp3d.o $(INC)/kelp_context.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_rte3d: $(SRC)/test_rte3d.f90 $(INC)/rte_sparse_matrices.o $(INC)/test_rte3d_mod.o $(INC)/mgmres.o $(INC)/rte3d.o $(INC)/test_kelp3d_mod.o $(INC)/prob.o $(INC)/fastgl.o $(INC)/sag.o $(INC)/utils.o $(INC)/kelp3d.o $(INC)/kelp_context.o $(INC)/hdf5_utils.o
	$(H5FC) $(BFLAGS) $^ -o $(BIN)/$@
### Old tests
old: test_interp test_rte2d test_vsf

test_interp: $(SRC)/test_interp.f90 $(INC)/utils.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_rte2d: $(SRC)/test_rte2d.f90 $(INC)/rte2d.o $(INC)/rte_core.o $(INC)/utils.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@
test_vsf: $(SRC)/test_vsf.f90 $(INC)/rte_core.o $(INC)/utils.o
	$(FC) $(BFLAGS) $^ -o $(BIN)/$@

################
# Object files #
################

$(INC)/sag.o: $(SRC)/sag.f90 $(INC)/utils.o $(INC)/fastgl.o
	$(FC) $(OFLAGS) $< -o $@
$(INC)/hdf5_utils.o: $(SRC)/hdf5_utils.f90 $(INC)/utils.o $(INC)/kelp_context.o
	$(H5FC) $(OFLAGS) $< -o $@
$(INC)/kelp_context.o: $(SRC)/kelp_context.f90 $(INC)/sag.o $(INC)/prob.o
	$(FC) $(OFLAGS) $< -o $@
$(INC)/light_context.o: $(SRC)/kelp_context.f90 $(INC)/kelp3d.o $(INC)/kelp_context.o $(INC)/sag.o $(INC)/prob.o
	$(FC) $(OFLAGS) $< -o $@
$(INC)/kelp3d.o: $(SRC)/kelp3d.f90 $(INC)/kelp_context.o
	$(FC) $(OFLAGS) $< -o $@
$(INC)/test_kelp3d_mod.o: $(SRC)/test_kelp3d_mod.f90 $(INC)/kelp3d.o $(INC)/hdf5_utils.o
	$(FC) $(OFLAGS) $< -o $@
$(INC)/test_rte3d_mod.o: $(SRC)/test_rte3d_mod.f90 $(INC)/test_kelp3d_mod.o $(INC)/rte3d.o $(INC)/kelp3d.o $(INC)/hdf5_utils.o 
	$(H5FC) $(OFLAGS) $< -o $@
$(INC)/rte_sparse_matrices.o: $(SRC)/rte_sparse_matrices.f90 $(INC)/sag.o $(INC)/kelp_context.o $(INC)/mgmres.o $(INC)/hdf5_utils.o
	$(H5FC) $(OFLAGS) $< -o $@
$(INC)/rte3d.o: $(SRC)/rte3d.f90 $(INC)/kelp_context.o $(INC)/rte_sparse_matrices.o
	$(FC) $(OFLAGS) $< -o $@
# Old

$(INC)/rte_core.o: $(SRC)/rte_core.f90 $(INC)/utils.o
	$(FC) $(OFLAGS) $< -o $@
$(INC)/utils.o: $(SRC)/utils.f90
	$(FC) $(OFLAGS) $< -o $@
$(INC)/rte2d.o: $(SRC)/rte2d.f90 $(INC)/rte_core.o
	$(FC) $(OFLAGS) $< -o $@
# External

$(INC)/fastgl.o: $(EXT)/fastgl.f90
	$(FC) $(OFLAGS) $< -o $@
$(INC)/mgmres.o: $(EXT)/mgmres.f90
	$(FC) $(OFLAGS) $< -o $@
$(INC)/prob.o: $(EXT)/prob.f90
	$(FC) $(OFLAGS) $< -o $@
# Utils

clean: rmo
	rm -f $(INC)/*.mod $(INC)/*.o $(BIN)/*

rmo: 
	rm -f $(BASE)/*.o $(BASE)/*.mod

ls:
	ls $(SRC) $(BIN) $(INC)
