#
# Process this file with autoconf to produce a configure script.
#
AC_INIT([GammaLib], [0.1pre], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR(src/GammaLib.hpp)

AM_INIT_AUTOMAKE([gnu no-dependencies])
AM_CONFIG_HEADER([config.h])
#AM_MAINTAINER_MODE

# things required by automake
# AC_ARG_PROGRAM
AC_PROG_MAKE_SET

# Check for which system.
# AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
# AC_LANG(C)
# AC_PROG_CC
# AC_PROG_CPP
AC_PROG_INSTALL
# AC_PROG_LN_S

# Library versioning (current:revision:age)
# See the libtool manual for an explanation
#
# gammalib-0.1     libgamma 0:0:0
GAMMALIB_LT_VERSION="0:0:0"
AC_SUBST(GAMMALIB_LT_VERSION)

# libtool test.
AC_PROG_LIBTOOL

# Check compiler features
# AC_TYPE_SIZE_T
# AC_C_CONST
# AC_C_VOLATILE
# AC_C_INLINE

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Check for range cheking option.
#RANGE_CHECK="no"
#AC_ARG_ENABLE(range_check, 
#             [  --enable-range-check    performs index range checking], 
#             [ RANGE_CHECK="$enableval"
#             ]) 
#if test "x$RANGE_CHECK" = "xyes"; then
#        AC_DEFINE([G_RANGE_CHECK], [1], 
#        [Define if range checking should be performed])
#fi

#
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
