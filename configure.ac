#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(ntsm, 1.2.1, cjustin@ds.dfci.harvard.edu)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

#see http://stackoverflow.com/questions/8215785/how-can-i-use-ac-revision-with-git
AC_REVISION([m4_esyscmd([./tools/configure.commit])])
AC_DEFINE([GIT_REVISION], "[m4_esyscmd([./tools/configure.commit])]", [git revision information])

AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Set compiler flags
AC_SUBST(CXXFLAGS,"-O3 -std=c++11 $CXXFLAGS")
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror $AM_CXXFLAGS")

#set dimension parameter (determined at compile time)
AC_ARG_VAR(NTSM_DIMS, dimensions to use with PCA)

# Find the absolute paths.
srcdir=$(cd $srcdir; pwd)

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Checks for typedefs, structures, and compiler characteristics.

AC_CONFIG_FILES([
		Makefile
		src/Makefile
])
AC_OUTPUT

