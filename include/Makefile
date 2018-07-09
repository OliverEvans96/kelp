# Kelp light model Makefile

# Project directories
BASE=..
INC=.
BIN=$(BASE)/bin
FTN=$(BASE)/code/fortran
SRC=$(FTN)/src
LIB=$(PWD)/../lib
EXT=$(FTN)/external
HDF5=$(FTN)/hdf5
TEST=$(FTN)/test
TEST_MOD=$(TEST)/mod
F2PY=$(FTN)/f2py
JTEST=$(TEST)/for_julia
PYDIR=$(BASE)/code/python
JULIA=$(BASE)/code/julia
# F2PY .o directory
F2PY_INC=$(PYDIR)/fortran_wrappers


## LIS (iterative solver)
# LIS location
LIS=$(BASE)/lis
LISFLAGS=-cpp -llis -lm -Wl,-rpath=$(LIB)
lisobj=$(LIB)/liblis.a

$(lisobj):
	git submodule update --init --progress # Get submodule if not present
	cd $(LIS) && ./configure --prefix=$(PWD)/.. --enable-shared --enable-omp --enable-fortran --enable-saamg --with-pic
	cd $(LIS) && make -j5
	cd $(LIS) && make install

lis: $(lisobj)

# Fortran Compiler
FC=ftn
# HDF5 Fortan Compiler
H5FC=h5fc
# MPI Fortran Compiler
MPIFC=mpif90

# Profiling (timing) flags
# Compile & run normally,
# Then call gprof `executable`
PFLAGS= -g -ffpe-trap=zero,invalid,overflow,underflow -ffpe-summary=all -fcheck=all -O0

# Optimize performance
#https://stackoverflow.com/questions/42386065/inlining-functions-in-fortran
# WARNING: Using -Ofast will break things due to allocating too large
# of arrays on the stack
# http://blog.purevirtual.net/2012/01/valgrind-client-switching-stacks.html
OPTFLAGS= -flto -fopenmp #-Ofast

# Doesn't work as shared? Executable gives SEGFAULT immediately.
# But that's okay, because this works.
FFLAGS=-fPIC $(LISFLAGS) -I$(INC) -L$(LIB) $(OPTFLAGS) $(PFLAGS) -Wall 
# Flags for F2PY
F2PYFLAGS=-fPIC $(OPTFLAGS) # $(PFLAGS) 
# Test flags
TESTMODFLAGS=-shared -fPIC $(FFLAGS)
# Export LD_LIBRARY_PATH for tests
TESTEXP=export LD_LIBRARY_PATH=$(INC):$(LD_LIBRARY_PATH)
TEST_FLAGS=-I/usr/include -L/usr/lib -lnetcdff

# All source files in each directory
f90_src = $(wildcard $(SRC)/*.f90)
ext_src = $(wildcard $(EXT)/*.f90)
hdf5_src = $(wildcard $(HDF5)/*.f90)
test_src = $(wildcard $(TEST)/*.f90)
test_mod_src = $(wildcard $(TEST_MOD)/*.f90)
jtest_src = $(wildcard $(JTEST)/*.f90)
f2py_src = $(wildcard $(F2PY)/*.f90)
# Names of .o files which they will produce
f90_inc = $(f90_src:$(SRC)/%.f90=%.o)
ext_inc = $(ext_src:$(EXT)/%.f90=%.o)
hdf5_inc = $(hdf5_src:$(HDF5)/%.f90=%.o)
test_mod_inc = $(test_mod_src:$(TEST_MOD)/%.f90=%.o)
jtest_inc = $(jtest_src:$(JTEST)/%.f90=%.so)
f2py_inc = $(f2py_src:$(F2PY)/%.f90=%.o)
# Other
test_bin = $(test_src:$(TEST)/%.f90=%)
f2py_targets = $(f2py_src:$(F2PY)/%.f90=%)

##########
# GROUPS #
##########

.PHONY: all lis clean cleanlis cleanall $(test_bin) $(f2py_targets)

# Default target makes all objects
all: $(f90_inc)

report:
#	@echo $(test_mod_src)
#	@echo $(test_mod_inc)
#	@echo $(test_bin)
	@echo $(f2py_targets)

# Make all f2py wrappers
f2py: $(f2py_targets)

#########
# RULES #
#########

# Declare which .o files depend on which .f90 files
# and provide the actual compilation command
$(f90_inc): %.o: $(SRC)/%.f90
	$(FC) -c $< $(FFLAGS) -o $@
$(ext_inc): %.o: $(EXT)/%.f90
	$(FC) -c $< $(FFLAGS) -o $@
$(hdf5_inc): %.o: $(HDF5)/%.f90
	$(H5FC) -c $< $(FFLAGS) -o $@
$(test_bin): %: $(TEST)/%.f90
	$(FC) $^ $(TEST_FLAGS) $(FFLAGS) -o $(BIN)/$@
$(test_mod_inc): %.o: $(TEST_MOD)/%.f90
	$(FC) -c $< $(FFLAGS) -o $@
$(jtest_inc): %.so: $(JTEST)/%.f90
	$(FC) $(TESTFLAGS) $< -o $@
# For actually importing in python
$(f2py_targets): %: $(F2PY)/%.f90
# TODO: This is sloppy & adds patchelf (Nix) dependency
	f2py -c $^ -m $@ --f90flags="$(F2PYFLAGS)" -llis -lm -L$(LIB)/../lib -I$(INC)
	patchelf --set-rpath $(LIB) $(PWD)/$@.*.so
	rm -f $(F2PY_INC)/$@.*.so
	ln -s $(PWD)/$@.*.so $(F2PY_INC)
# To test outside of python directly in fortran
$(f2py_inc): %.o: $(F2PY)/%.f90
	$(FC) -c $< $(FFLAGS) -o $@


################
# Object files #
################

# FC
rte_core.o: utils.o
sag.o: utils.o fastgl.o
kelp_context.o: sag.o prob.o
kelp3d.o: kelp_context.o
rte_sparse_matrices.o: sag.o kelp_context.o mgmres.o $(lisobj)
light_context.o: sag.o rte_sparse_matrices.o
rte3d.o: kelp_context.o rte_sparse_matrices.o light_context.o
asymptotics.o: kelp_context.o rte_sparse_matrices.o light_context.o
#pykelp3d_wrap.o: kelp3d.o
#pyrte3d_wrap.o: rte3d.o
#pyasymptotics_wrap.o: asymptotics.o
light_interface.o: asymptotics.o rte3d.o kelp3d.o

# H5FC
hdf5_utils.o: utils.o kelp_context.o

# Test mods
test_kelp3d_mod.o: kelp3d.o #hdf5_utils.o
test_rte3d_mod.o: test_kelp3d_mod.o rte3d.o kelp3d.o light_context.o #hdf5_utils.o
test_grid_mod.o: sag.o
test_asymptotics.o: asymptotics.o light_context.o rte_sparse_matrices.o mgmres.o rte3d.o
test_asymptotics_mod.o: asymptotics.o prob.o fastgl.o sag.o utils.o rte3d.o kelp_context.o light_context.o rte_sparse_matrices.o asymptotics.o
# F2PY test w/o python
pykelp3d_wrap.o: light_interface.o asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o

# F2PY Test
test_pykelp3d_wrap: pykelp3d_wrap.o asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o
nc_pykelp3d_wrap: pykelp3d_wrap.o asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o
test_pyasymptotics_wrap: asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o

########
# F2PY #
########

pykelp3d_wrap: light_interface.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o mgmres.o rte_sparse_matrices.o light_context.o rte3d.o asymptotics.o

#########
# Tests #
#########

#tests: test_context test_gl test_asymptotics test_grid test_gmres #test_prob test_kelp_3d

tests: $(test_bin)


# Julia Tests
test: $(jtest_inc)
	$(TESTEXP); cd $(JULIA); julia tests.jl

# TESTS for use with Julia
test_grid.so: fastgl.o sag.o kelp_context.o prob.o utils.o
test_asymptotics.so: asymptotics.o light_context.o rte_sparse_matrices.o mgmres.o rte3d.o kelp_context.o fastgl.o prob.o sag.o utils.o

test_test_asymptotics: test_asymptotics.so

test_traverse: test_asymptotics.o asymptotics.o light_context.o rte_sparse_matrices.o mgmres.o rte3d.o kelp_context.o fastgl.o prob.o sag.o utils.o
# Old Fortran tests
test_context: prob.o utils.o kelp_context.o
test_gl: fastgl.o
test_grid: test_grid_mod.o kelp_context.o fastgl.o prob.o sag.o utils.o
test_asymptotics: test_asymptotics_mod.o asymptotics.o light_context.o rte_sparse_matrices.o mgmres.o rte3d.o kelp_context.o fastgl.o prob.o sag.o utils.o
test_gmres: mgmres.o utils.o hdf5_utils.o kelp_context.o sag.o fastgl.o prob.o
test_prob: prob.o prob.o
test_kelp3d: test_kelp3d_mod.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o #hdf5_utils.o
test_rte3d: test_rte3d_mod.o #rte_sparse_matrices.o test_rte3d_mod.o mgmres.o rte3d.o test_kelp3d_mod.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o #hdf5_utils.o
test_pyrte3d_wrap: pyrte3d_wrap.o rte_sparse_matrices.o mgmres.o rte3d.o test_kelp3d_mod.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o #hdf5_utils.o
#test_pyasymptotics_wrap: pykelp3d_wrap #asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o test_kelp3d_mod.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o #hdf5_utils.o
test_angular_integral: light_interface.o pyasymptotics_wrap.o asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o test_kelp3d_mod.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o #hdf5_utils.o
test_light_interface: light_interface.o asymptotics.o rte_sparse_matrices.o mgmres.o rte3d.o prob.o fastgl.o sag.o utils.o kelp3d.o kelp_context.o light_context.o #hdf5_utils.o

test_interp: utils.o
test_rte2d: rte2d.o rte_core.o utils.o
test_vsf: rte_core.o utils.o

#############
# UTILITIES #
#############

clean:
	rm -f *.mod *.so *.o $(F2PY_INC)/*.so $(BIN)/*

cleanlis:
	cd $(LIS) && make clean
	rm -f $(LIB)/*.a $(LIB)/*.la $(LIB)/*.so
	rm -f $(INC)/*.h

cleanall: clean cleanlis
