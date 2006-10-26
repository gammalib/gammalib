# Process this file with autoconf to produce a configure script.

AC_INIT([gamma], [0.1], [knodlseder@cesr.fr])
# AC_CONFIG_SRCDIR(src/gammalib.hpp)
AC_CONFIG_SRCDIR(src/GException.hpp)

AM_INIT_AUTOMAKE([gnu no-dependencies])
AM_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE

# things required by automake
# AC_ARG_PROGRAM
AC_PROG_MAKE_SET

# Check for which system.
# AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)
# AC_LANG(C)
# AC_PROG_CC
# AC_PROG_CPP
# AC_PROG_INSTALL
# AC_PROG_LN_S

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
# AC_TYPE_SIZE_T
# AC_C_CONST
# AC_C_VOLATILE
# AC_C_INLINE


AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
