# Kelp light model Makefile

# Project directories
BASE=..
BIN=$(BASE)/bin
FTN=$(BASE)/code/fortran
SRC=$(FTN)/src
EXT=$(FTN)/external
TEST=$(FTN)/test
F2PY=$(FTN)/f2py
JTEST=$(TEST)/for_julia
INC=.
PYDIR=$(BASE)/code/python
JULIA=$(BASE)/code/julia
# F2PY .so / .py directory

# Fortran Compiler
FC=gfortran
# HDF5 Fortan Compiler
H5FC=h5fc

# Profiling (timing) flags
# Add to OFLAGS and BFLAGS
# (Should replace CFLAGS in OFLAGS)
# Compile & run normally,
# Then call gprof `executable`
PFLAGS= -g -O0 -ffpe-trap=zero,invalid,overflow,underflow -ffpe-summary=all #-fcheck=all

# Optimize performance
#https://stackoverflow.com/questions/42386065/inlining-functions-in-fortran
#OPTFLAGS=-Ofast -flto

# Doesn't work as shared? Executable gives SEGFAULT immediately.
# But that's okay, because this works.
FFLAGS = $(PFLAGS) -fPIC #-fopenmp
# Flags for F2PY
F2PYFLAGS= -fPIC # $(PFLAGS) $(OPTFLAGS) #$(OMPFLAGS)
# Test flags
TESTFLAGS= -shared -fPIC $(FFLAGS)
# Export LD_LIBRARY_PATH for tests
TESTEXP=export LD_LIBRARY_PATH=$(INC):$(LD_LIBRARY_PATH)

# All source files in each directory
f90_src = $(wildcard $(SRC)/*.f90)
ext_src = $(wildcard $(EXT)/*.f90)
test_src = $(wildcard $(TEST)/*.f90)
jtest_src = $(wildcard $(JTEST)/*.f90)
f2py_src = $(wildcard $(F2PY)/*.f90)
# Names of .so files which they will produce
f90_inc = $(f90_src:$(SRC)/%.f90=%.so)
ext_inc = $(ext_src:$(EXT)/%.f90=%.so)
test_bin = $(test_src:$(TEST)/%.f90=%)
jtest_inc = $(jtest_src:$(JTEST)/%.f90=%.so)
f2py_targets = $(f2py_src:$(F2PY)/%.f90=%)
.PHONY: $(test_bin) $(f2py_targets)

# Default target makes all objects
all: $(f90_inc)

# Declare which .so files depend on which .f90 files
# and provide the actual compilation command
$(f90_inc): %.so: $(SRC)/%.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(ext_inc): %.so: $(EXT)/%.f90
	$(FC) -c $(FFLAGS) $< -o $@
$(test_bin): %: $(TEST)/%.f90
	$(FC) $(FFLAGS) $^ -o $(BIN)/$@
$(jtest_inc): %.so: $(JTEST)/%.f90
	$(FC) $(TESTFLAGS) $^ -o $@
$(f2py_targets): %: $(F2PY)/%.f90
	f2py -m $@ -c $^

################
# Object files #
################

# FC
rte_core.so: utils.so
sag.so: utils.so fastgl.so
kelp_context.so: sag.so prob.so
kelp3d.so: kelp_context.so
rte_sparse_matrices.so: sag.so kelp_context.so mgmres.so
light_context.so: sag.so rte_sparse_matrices.so
rte3d.so: kelp_context.so rte_sparse_matrices.so light_context.so
asymptotics.so: kelp_context.so rte_sparse_matrices.so light_context.so
#pykelp3d_wrap.so: kelp3d.so
#pyrte3d_wrap.so: rte3d.so
#pyasymptotics_wrap.so: asymptotics.so
light_interface.so: asymptotics.so rte3d.so kelp3d.so

# Test mods
test_kelp3d_mod.so: kelp3d.so #hdf5_utils.so
test_rte3d_mod.so: test_kelp3d_mod.so rte3d.so kelp3d.so light_context.so #hdf5_utils.so
test_grid_mod.so: sag.so
test_asymptotics.so: asymptotics.so light_context.so rte_sparse_matrices.so mgmres.so rte3d.so
test_asymptotics_mod.so: asymptotics.so prob.so fastgl.so sag.so utils.so rte3d.so kelp_context.so light_context.so rte_sparse_matrices.so asymptotics.so

# H5FC
hdf5_utils.so: utils.so kelp_context.so
	$(H5FC) -c $(OFLAGS) $< -o $@

########
# F2PY #
########

pykelp3d_wrap: pykelp3d_wrap.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so

pyrte3d_wrap: pyrte3d_wrap.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so mgmres.so rte_sparse_matrices.so light_context.so rte3d.so

pyasymptotics_wrap: prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so mgmres.so rte_sparse_matrices.so light_context.so rte3d.so asymptotics.so

#########
# Tests #
#########

#tests: test_context test_gl test_asymptotics test_grid test_gmres #test_prob test_kelp_3d

tests: $(test_bin)


# Julia Tests
test: $(jtest_inc)
	$(TESTEXP); julia $(JULIA)/tests.jl

# TESTS for use with Julia
test_grid.so: fastgl.so sag.so kelp_context.so prob.so utils.so
test_asymptotics.so: asymptotics.so light_context.so rte_sparse_matrices.so mgmres.so rte3d.so kelp_context.so fastgl.so prob.so sag.so utils.so

test_traverse: test_asymptotics.so asymptotics.so light_context.so rte_sparse_matrices.so mgmres.so rte3d.so kelp_context.so fastgl.so prob.so sag.so utils.so
# Old Fortran tests
test_context: prob.so utils.so kelp_context.so
test_gl: fastgl.so
test_grid: test_grid_mod.so kelp_context.so fastgl.so prob.so sag.so utils.so
test_asymptotics: test_asymptotics_mod.so asymptotics.so light_context.so rte_sparse_matrices.so mgmres.so rte3d.so kelp_context.so fastgl.so prob.so sag.so utils.so
test_gmres: mgmres.so utils.so hdf5_utils.so kelp_context.so sag.so fastgl.so prob.so
test_prob: prob.so prob.so
test_kelp3d: test_kelp3d_mod.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so #hdf5_utils.so
test_rte3d: rte_sparse_matrices.so test_rte3d_mod.so mgmres.so rte3d.so test_kelp3d_mod.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so light_context.so #hdf5_utils.so
test_pyrte3d_wrap: pyrte3d_wrap.so rte_sparse_matrices.so mgmres.so rte3d.so test_kelp3d_mod.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so light_context.so #hdf5_utils.so
test_pyasymptotics_wrap: pyasymptotics_wrap.so asymptotics.so rte_sparse_matrices.so mgmres.so rte3d.so test_kelp3d_mod.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so light_context.so #hdf5_utils.so
test_angular_integral: light_interface.so pyasymptotics_wrap.so asymptotics.so rte_sparse_matrices.so mgmres.so rte3d.so test_kelp3d_mod.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so light_context.so #hdf5_utils.so
test_light_interface: light_interface.so asymptotics.so rte_sparse_matrices.so mgmres.so rte3d.so prob.so fastgl.so sag.so utils.so kelp3d.so kelp_context.so light_context.so #hdf5_utils.so

test_interp: utils.so
test_rte2d: rte2d.so rte_core.so utils.so
test_vsf: rte_core.so utils.so

#############
# UTILITIES #
#############

.PHONY: clean

clean:
	rm -f *.mod *.so *.py $(BIN)/*
