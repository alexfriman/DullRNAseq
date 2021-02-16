#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([DullRNAseq], [0.01], [alex@friman.online])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/genomeparser.cpp])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
])
# Checks for programs. 
AC_PROG_CXX([mpic++ mpic++.openmpi mpicxx mpicxx.openmpi mpig++])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.

XXX_PROGRAM="
#include <fstream>
#include <cstring>
#include <stdio.h>
#include <string>
#include <vector> 
#include <mpi.h>
int main(void) {
    int x=0;
    return x;
}
"

AC_LANG([C++])
AC_CHECK_FUNCS([memset])
AM_INIT_AUTOMAKE([foreign subdir-objects])
#AM_INIT_AUTOMAKE([foreign])
#AUTOMAKE_OPTIONS = foreign
AC_MSG_CHECKING([whether $CXX finds fstream cstring string stdio.h vector mpi.h])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([$XXX_PROGRAM])],[AC_MSG_RESULT([yes])],[AC_MSG_ERROR([Failed to find header files])])


AC_OUTPUT
