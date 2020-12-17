bin_PROGRAMS = heatsolver
AM_CXXFLAGS = $(MASA_CXXFLAGS) $(GRVY_CFLAGS) $(HDF5_CFLAGS)

if PETSC_ENABLED
   include $(PETSC_DIR)/lib/petsc/conf/variables
   AM_CXXFLAGS += -DINCLUDE_PETSC $(PETSC_CC_INCLUDES)
endif

LIBS = $(MASA_LIBS) $(GRVY_LIBS) $(HDF5_LIBS)
heatsolver_SOURCES = heat.cc masa_arrays.cc masa_arrays.h matrices.cc matrices.h numerical_solve.cc numerical_solve.h error_norm.cc error_norm.h parameter_check.cc parameter_check.h node_val_check.cc node_val_check.h parser.cc parser.h hd_prep.cc hd_prep.h hdf5_output.cc hdf5_output.h petsc_error.cc petsc_error.h
if PETSC_ENABLED
LIBS += $(PETSC_LIB)
heatsolver_SOURCES += petsc_setup.cc petsc_setup.h petsc_masa.cc petsc_masa.h petsc_solve.cc petsc_solve.h
endif
