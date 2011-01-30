# Configuration for GammaLib - a versatile toolbox for high-level analysis
#                              of astronomical gamma-ray data
#
# Copyright (C) 2006-2011
#
# GammaLib is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3, or (at your option) any
# later version.
#
# Process this file with autoconf to produce a configure script

AC_INIT([GammaLib], [0.4.4], [knodlseder@cesr.fr], [gammalib], [http://gammalib.sourceforge.net/])
AC_CONFIG_SRCDIR([src/support/GTools.cpp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu no-dependencies])
AC_CONFIG_HEADERS([config.h])


dnl =========================================
dnl Library versioning (current:revision:age)
dnl See the libtool manual for an explanation
dnl
dnl gammalib-0.1     libgamma 0:0:0
dnl gammalib-0.2     libgamma 0:0:0
dnl gammalib-0.3     libgamma 0:0:0
dnl gammalib-0.4     libgamma 0:0:0
dnl =========================================
GAMMALIB_LT_VERSION="0:0:0"
AC_SUBST(GAMMALIB_LT_VERSION)


dnl ==========
dnl Set prefix
dnl ==========
AC_PREFIX_DEFAULT([/usr/local/gamma])
if test "x$prefix" = xNONE; then
  AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], ["/usr/local/gamma"],
                     [Location where library is installed])
else
  AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${prefix}",
                     [Location where library is installed])
fi


dnl ===========================
dnl things required by automake
dnl ===========================
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET


dnl =========================
dnl Check for system services
dnl =========================
AC_CANONICAL_HOST


dnl =======================
dnl Checks for header files
dnl =======================
AC_STDC_HEADERS


dnl =============================================================
dnl Checks for typedefs, structures, and compiler characteristics
dnl =============================================================
AC_TYPE_SIZE_T
AC_C_CONST
AC_C_INLINE


dnl ===================
dnl Checks for programs
dnl ===================
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL


dnl =================================
dnl Checks for debugging or profiling
dnl =================================
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE([debug],
              AS_HELP_STRING([--enable-debug],
                             [turn on debugging (default=no)]),
              [enable_debug="$enableval"],
              [enable_debug="no"])
if test "x$enable_debug" = "xyes"; then
  CFLAGS="$CFLAGS -g"
  CXXFLAGS="$CXXFLAGS -g"
  FFLAGS="$FFLAGS -g"
  AC_DEFINE([G_DEBUG], [1], [Define if debugging is turned on])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi
AC_MSG_CHECKING(whether to enable profiling)
AC_ARG_ENABLE([profiling],
              AS_HELP_STRING([--enable-profiling],
                             [turn on profiling (default=no)]),
              [enable_profiling="$enableval"],
              [enable_profiling="no"])
if test "x$enable_profiling" = "xyes"; then
  CFLAGS="$CFLAGS -pg"
  CXXFLAGS="$CXXFLAGS -pg"
  FFLAGS="$FFLAGS -pg"
  AC_DEFINE([G_PROFILE], [1], [Define if profiling is turned on])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi


dnl ============================
dnl Checks for library functions
dnl ============================


dnl =====================================
dnl Range checking option (G_RANGE_CHECK)
dnl =====================================
AC_MSG_CHECKING(whether to enable range checking)
AC_ARG_ENABLE([range_check],
              AS_HELP_STRING([--enable-range-check],
                             [performs index range checking (default=yes)]),
             [RANGE_CHECK="$enableval"],
             [RANGE_CHECK="yes"])
if test "x$RANGE_CHECK" = "xyes"; then
  AC_DEFINE([G_RANGE_CHECK], [1], [Define if range checking should be performed])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi


dnl ====================================
dnl Small memory option (G_SMALL_MEMORY)
dnl ====================================
AC_MSG_CHECKING(whether to enable small memory option)
AC_ARG_ENABLE([small_memory],
              AS_HELP_STRING([--enable-small-memory],
                             [optimizes code for small memory (default=yes)]),
             [SMALL_MEMORY="$enableval"],
             [SMALL_MEMORY="yes"])
if test "x$SMALL_MEMORY" = "xyes"; then
  AC_DEFINE([G_SMALL_MEMORY], [1], [Define if small memory optimization is selected])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi


dnl ========================================================================
dnl Checks for python binding. Python binding requires Python and either the
dnl files gammalib_wrap.cpp and gammalib.py or SWIG to create these files.
dnl ========================================================================
# Checks for Python (HAVE_PYTHON)
AM_PATH_PYTHON([2.2],
               [],
               [AC_MSG_WARN([Python not found. Python is required to build python binding. Python can be obtained from http://www.python.org])])

# Python binding enabled?
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding (default=yes)]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    # Checks for gammalib_wrap.cpp and gammalib.py
    AC_CHECK_FILES([pyext/gammalib_wrap.cpp pyext/gammalib.py],
                   [ac_has_python_wrap="yes"],
                   [ac_has_python_wrap="no"])
    if test "x$ac_has_python_wrap" = "xyes"; then
      AC_MSG_NOTICE([Build Python binding using shipped gammalib_wrap.cpp and gammalib.py files])
    else
      # Checks for SWIG (HAVE_SWIG)
      AC_PATH_PROG([SWIG], [swig])
      if test ! "$SWIG"; then
        AC_MSG_WARN([Neither pyext/gammalib_wrap.cpp and pyext/gammalib.py nor swig has not been found])
        AC_MSG_WARN([Python binding for gammalib cannot be built])
        AC_MSG_WARN([Either use a gammalib distribution with gammalib_wrap.cpp and gammalib.py shipped or install swig on your system])
        AC_MSG_WARN([swig can be obtained from http://www.swig.org/])
      else
        AC_MSG_NOTICE([Build Python binding using swig])
      fi
    fi
    # Set host specific build prefix (needed for Solaris bug)
    PYTHON_BUILD_PREFIX=
    case $host in
      sparc*-sun-solaris*)  PYTHON_BUILD_PREFIX="CC=CC" ;;
      *)                    PYTHON_BUILD_PREFIX="" ;;
    esac                    
    AC_SUBST(PYTHON_BUILD_PREFIX)
  else
    AC_MSG_WARN([Python binding for gammalib cannot be built as Python seems not to be installed on your system])
  fi
else
  AC_MSG_NOTICE([Python binding for gammalib will not be built])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi
AM_CONDITIONAL(ENABLE_PYTHON_BINDING, test "x$ac_enable_python_binding" = "xyes")
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")
AM_CONDITIONAL(HAS_PYTHON_WRAP, test "x$ac_has_python_wrap" = "xyes")
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")


dnl ==========================
dnl Checks for cfitsio library
dnl ==========================
# Set host specific cfitsio dependencies
cfitsio_deps=
case $host in
  sparc*-sun-solaris*)  cfitsio_deps="-lm -lsocket" ;;
  *)                    cfitsio_deps="-lm" ;;
esac                    

# Check if cfitsio is available
AC_CHECK_LIB([cfitsio], [ffpss],
             [LIBS="${LIBS} -lcfitsio ${cfitsio_deps}"
              AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])],
             [], [${cfitsio_deps}])
AM_CONDITIONAL(WITH_CFITSIO, test "x$ac_cv_lib_cfitsio_ffpss" = "xyes")

# If cfitsio is available then get library path for unit testing
if test "x$ac_cv_lib_cfitsio_ffpss" = "xyes"; then
  for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
    test -f $i/libcfitsio.a -o -f $i/libcfitsio.so -o -f $i/libcfitsio.dylib && LIBDIR_CFITSIO=$i
  done
  AC_SUBST(LIBDIR_CFITSIO)
fi

# Check if we want to use cfitsio
#AC_ARG_WITH([cfitsio],
#            [AS_HELP_STRING([--with-cfitsio],
#                            [Use CFITSIO library (needed to support FITS files)])],
#            [],
#            [with_cfitsio=yes])
#AC_ARG_WITH([cfitsio-lib-dir],
#            AC_HELP_STRING([--with-cfitsio-lib-dir],
#                           [CFITSIO library directory [[/usr/(local/)lib(64)]]]),
#            [cfitsio_lib_dir=$withval])
#
# If cfitsio is requested then add cfitsio library path to LDFLAGS
#if test "x$with_cfitsio" = "xyes"; then
#  if test "x$cfitsio_lib_dir" != "x"; then
#    test -f $cfitsio_lib_dir/libcfitsio.a -o -f $cfitsio_lib_dir/libcfitsio.so -o -f $cfitsio_lib_dir/libcfitsio.dylib && LIBDIR_CFITSIO=$cfitsio_lib_dir
#  else
#    for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
#      test -f $i/libcfitsio.a -o -f $i/libcfitsio.so -o -f $i/libcfitsio.dylib && LIBDIR_CFITSIO=$i
#    done
#    for i in /usr/local/include /opt/local/include /usr/include ; do
#      test -f $i/fitsio.h && INCDIR_CFITSIO=$i
#    done
#  fi
#  if test -z "$LIBDIR_CFITSIO"; then
#    if test "x$cfitsio_lib_dir" != "x";then
#      AC_MSG_ERROR(Cannot find CFITSIO library under $cfitsio_lib_dir)
#    else
#      AC_MSG_ERROR(Cannot find CFITSIO library. Use --with-cfitsio-lib-dir= to specify non-default path or --without-cfitsio to compile with CFITSIO library.)
#    fi
#  fi
#  LDFLAGS="$LDFLAGS -L$LIBDIR_CFITSIO"
#  CPPFLAGS="$CPPFLAGS -I$LIBDIR_CFITSIO/../include"
#  AC_SUBST(LIBDIR_CFITSIO)
#fi
#
# Check for cfitsio library
#if test "x$with_cfitsio" = "xyes"; then
#  AC_CHECK_LIB([cfitsio], [ffpss],
#               [LIBS="${LIBS} -lcfitsio -lm"
#                AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])],
#               [AC_MSG_FAILURE([--with-cfitsio was given, but test for cfitsio failed])],
#               [-lm])
#fi
#AM_CONDITIONAL(WITH_CFITSIO, test "x$with_cfitsio" = "xyes")


dnl ===========================
dnl Checks for readline library
dnl ===========================

# Set host specific readline dependencies
readline_deps=
case $host in
  sparc*-sun-solaris*)  readline_deps="-lncurses" ;;
  *)                    readline_deps="-lncurses" ;;
esac                    

# Check if readline is available
AC_CHECK_LIB([readline], [main],
             [LIBS="${LIBS} -lreadline ${readline_deps}"
              AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])],
             [], [${readline_deps}])
AM_CONDITIONAL(WITH_READLINE, test "x$ac_cv_lib_readline_main" = "xyes")

# If readline is available then get library path for unit testing
if test "x$ac_cv_lib_readline_main" = "xyes"; then
  for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
    test -f $i/libreadline.a -o -f $i/libreadline.so -o -f $i/libreadline.dylib && LIBDIR_READLINE=$i
  done
  AC_SUBST(LIBDIR_READLINE)
fi

# Checks for readline
#AC_ARG_WITH([readline],
#            [AS_HELP_STRING([--with-readline],
#            [support fancy command line editing @<:@default=check@:>@])],
#            [],
#            [with_readline=yes])
#AC_ARG_WITH([readline-lib-dir],
#            AC_HELP_STRING([--with-readline-lib-dir],
#                           [readline library directory [[/usr/(local/)lib(64)]]]),
#            [readline_lib_dir=$withval])
#
# If readline is requested then add readline library path to LDFLAGS
#if test "x$with_readline" = "xyes"; then
#  if test "x$readline_lib_dir" != "x"; then
#    test -f $readline_lib_dir/libreadline.a -o -f $readline_lib_dir/libreadline.so -o -f $readline_lib_dir/libreadline.dylib && LIBDIR_READLINE=$readline_lib_dir
#  else
#    for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
#      test -f $i/libreadline.a -o -f $i/libreadline.so -o -f $i/libreadline.dylib && LIBDIR_READLINE=$i
#    done
#  fi
#  if test -z "$LIBDIR_READLINE"; then
#    if test "x$readline_lib_dir" != "x";then
#      AC_MSG_ERROR(Cannot find readline library under $readline_lib_dir)
#    else
#      AC_MSG_ERROR(Cannot find readline library. Use --with-readline-lib-dir= to specify non-default path or --without-readline to compile with readline library.)
#    fi
#  fi
#  LDFLAGS="$LDFLAGS -L$LIBDIR_READLINE"
#  AC_SUBST(LIBDIR_READLINE)
#fi
#
# Check for readline library
#if test "x$with_readline" = "xyes"; then
#  AC_CHECK_LIB([readline], [main],
#               [LIBS="${LIBS} -lreadline -lncurses"
#                AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])],
#               [AC_MSG_FAILURE([--with-readline was given, but test for readline failed])],
#               [-lncurses])
#fi
#AM_CONDITIONAL(WITH_READLINE, test "x$with_readline" = "xyes")


dnl =========================================
dnl Checks for instrument specific interfaces
dnl =========================================
# Checks for the existence of instrument specific interfaces
AC_CHECK_FILE([inst/Makefile.am],
              [with_inst=yes],
              [with_inst=no])
AM_CONDITIONAL([WITH_INST], [test "x$with_inst" = "xyes"])

# Checks for MWL interface
AC_ARG_WITH([mwl],
            [AS_HELP_STRING([--with-mwl],
                            [compile in multi-wavelength interface [default=yes]])],
            [],
            [with_mwl=yes])
if test "x$with_mwl" = "xyes"; then
  AC_CHECK_FILE([inst/mwl/Makefile.am],
                [with_mwl=yes],
                [with_mwl=no])
fi
AM_CONDITIONAL(WITH_INST_MWL, test "x$with_mwl" = "xyes")

# Checks for CTA interface
AC_ARG_WITH([cta],
            [AS_HELP_STRING([--with-cta],
                            [compile in CTA specific interface [default=yes]])],
            [],
            [with_cta=yes])
if test "x$with_cta" = "xyes"; then
  AC_CHECK_FILE([inst/cta/Makefile.am],
                [with_cta=yes],
                [with_cta=no])
fi
AM_CONDITIONAL(WITH_INST_CTA, test "x$with_cta" = "xyes")

# Checks for LAT interface
AC_ARG_WITH([lat],
            [AS_HELP_STRING([--with-lat],
                            [compile in LAT specific interface [default=yes]])],
            [],
            [with_lat=yes])
if test "x$with_lat" = "xyes"; then
  AC_CHECK_FILE([inst/lat/Makefile.am],
                [with_lat=yes],
                [with_lat=no])
fi
AM_CONDITIONAL(WITH_INST_LAT, test "x$with_lat" = "xyes")


dnl ======================================
dnl Set list of output files to be created
dnl ======================================
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/support/Makefile
                 src/linalg/Makefile
                 src/numerics/Makefile
                 src/fits/Makefile
                 src/xml/Makefile
                 src/sky/Makefile
                 src/opt/Makefile
                 src/obs/Makefile
                 src/model/Makefile
                 src/app/Makefile
                 include/Makefile
                 test/Makefile
                 pyext/Makefile
                 pyext/setup.py
                 gammalib.pc])
if test "x$with_inst" = "xyes"; then
  AC_CONFIG_FILES([inst/Makefile])
fi
if test "x$with_mwl" = "xyes"; then
  AC_CONFIG_FILES([inst/mwl/Makefile])
fi
if test "x$with_cta" = "xyes"; then
  AC_CONFIG_FILES([inst/cta/Makefile])
fi
if test "x$with_lat" = "xyes"; then
  AC_CONFIG_FILES([inst/lat/Makefile])
fi
AC_OUTPUT


dnl ===========================
dnl Print configuration summary
dnl ===========================
echo
echo "  GammaLib configuration summary"
echo "  =============================="
if test "x$ac_cv_lib_cfitsio_ffpss" = "xyes"; then
  echo "  * FITS I/O support             (yes)"
else
  echo "  - FITS I/O support             (no)    no cfitsio library found"
fi
if test "x$ac_cv_lib_readline_main" = "xyes"; then
  echo "  * Readline support             (yes)"
else
  echo "  - Readline support             (no)    no readline library found"
fi
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    if test "x$ac_has_python_wrap" = "xyes"; then
      echo "  * Python binding               (yes)"
    else
      if test ! "$SWIG"; then
        echo "  - Python binding               (no)    swig not installed"
      else
        echo "  * Python binding               (yes)   created by swig"
      fi
    fi
  else
    echo "  - Python binding               (no)    python not installed"
  fi
else
  echo "  - Python binding               (no)"
fi
if test "x$with_mwl" = "xyes"; then
  echo "  * Multiwavelength interface    (yes)"
else
  echo "  - Multiwavelength interface    (no)"
fi
if test "x$with_lat" = "xyes"; then
  echo "  * Fermi-LAT interface          (yes)"
else
  echo "  - Fermi-LAT interface          (no)"
fi
if test "x$with_cta" = "xyes"; then
  echo "  * CTA interface                (yes)"
else
  echo "  - CTA interface                (no)"
fi
echo
echo "Now type 'make'"
echo
