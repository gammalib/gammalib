dnl =============================================================
dnl Process this file with autoconf to produce a configure script
dnl =============================================================
AC_INIT([GammaLib], [0.2], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR(src/Makefile.am)

AM_INIT_AUTOMAKE([gnu no-dependencies])
AM_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE


dnl ===========================
dnl things required by automake
dnl ===========================
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET


dnl ======================
dnl Check for which system
dnl ======================
dnl AC_CANONICAL_HOST


dnl ===================
dnl Checks for programs
dnl ===================
AC_PROG_CXX
dnl AC_LANG(C)
dnl AC_PROG_CC
dnl AC_PROG_CPP
AC_PROG_INSTALL
dnl AC_PROG_LN_S


dnl ========================================= 
dnl Library versioning (current:revision:age)
dnl See the libtool manual for an explanation
dnl
dnl gammalib-0.1     libgamma 0:0:0
dnl gammalib-0.2     libgamma 0:0:0
dnl =========================================
GAMMALIB_LT_VERSION="0:0:0"
AC_SUBST(GAMMALIB_LT_VERSION)


dnl ============
dnl libtool test
dnl ============
AC_PROG_LIBTOOL


dnl =======================
dnl Check compiler features
dnl =======================
dnl AC_TYPE_SIZE_T
dnl AC_C_CONST
dnl AC_C_VOLATILE
dnl AC_C_INLINE


dnl =======================
dnl Checks for header files
dnl =======================


dnl =============================================================
dnl Checks for typedefs, structures, and compiler characteristics
dnl =============================================================


dnl ============================
dnl Checks for library functions
dnl ============================


dnl ==============================
dnl Check for range cheking option
dnl ==============================
RANGE_CHECK="no"
AC_ARG_ENABLE(range_check, 
             [  --enable-range-check    performs index range checking], 
             [ RANGE_CHECK="$enableval" ]) 
if test "x$RANGE_CHECK" = "xyes"; then
        AC_DEFINE([G_RANGE_CHECK], [1], [Define if range checking should be performed])
fi


dnl ===============
dnl Checks for SWIG
dnl ===============
AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")


dnl =================
dnl Checks for Python
dnl =================
AM_PATH_PYTHON([2.2],
               [],
               [AC_MSG_WARN([Python not found. Python is required to build soothsayer python binding. Python can be obtained from http://www.python.org])])
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")

AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding (default=yes]),
                             [ac_enable_python_binding=no],
                             [ac_enable_python_binding=yes])
if test "x$ac_enable_python_binding" = "xyes"
then
    if test ! "$SWIG" -a ! "$PYTHON"
    then
        AC_MSG_WARN([Python binding for gammalib cannot be built. Ensure that SWIG and Python packages are available.])
    fi
else
    AC_MSG_NOTICE([Python binding for gammalib will not be built.])
    AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi
AM_CONDITIONAL(ENABLE_PYTHON_BINDING, test "x$ac_enable_python_binding" = "xyes")


dnl ==================
dnl Checks for cfitsio
dnl ==================
#AC_SEARCH_LIBS(ffopen, cfitsio,
AC_SEARCH_LIBS(ffpss, cfitsio,
    AC_DEFINE([G_CFITSIO], [1], [Define if cfitsio is available]), ,[-lm])


dnl ===
dnl
dnl ===
AC_OUTPUT([Makefile src/Makefile include/Makefile test/Makefile pyext/Makefile pyext/setup.py])
