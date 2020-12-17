#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([heat-equation], [1.0], [aleccarruthers@utexas.edu])
AC_CONFIG_SRCDIR([src/heat.cc])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
#AC_PROG_CC

# Check for code coverage
AX_CODE_COVERAGE

# Checks for libraries.
AX_PATH_MASA(0.30,[yes],AC_MSG_ERROR([Masa cannot be found or version number does not meet the requirement]))
AX_PATH_HDF5(1.8.0,[yes],AC_MSG_ERROR([hdf5 cannot be found or version number does not meet the requirement]))
AX_PATH_GRVY(0.32,[yes],AC_MSG_ERROR([Grvy cannot be found or version number does not meet the requirement]))
AX_PATH_PETSC([no])
# Checks for header files.
#AC_CHECK_HEADERS([string.h])
AM_INIT_AUTOMAKE


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
#AC_CHECK_FUNCS([pow sqrt])

AC_OUTPUT(Makefile src/Makefile tests/Makefile)
