# Configuration for GammaLib - a versatile toolbox for high-level analysis
#                              of astronomical gamma-ray data
#
# Copyright (C) 2006-2011
#
# GammaLib is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any
# later version.
#
# Process this file with autoconf to produce a configure script

AC_INIT([GammaLib], [0.4.6], [knodlseder@cesr.fr], [gammalib], [http://gammalib.sourceforge.net/])
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
  gammalib_prefix="/usr/local/gamma"
else
  gammalib_prefix=$prefix
fi
AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${gammalib_prefix}", [Location where library is installed])


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
dnl AC_TYPE_LONG_LONG_INT
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
dnl The following conditionals are set:
dnl  HAVE_PYTHON           - Python >= 2.2 is available
dnl  HAVE_SWIG             - swig is available
dnl  ENABLE_PYTHON_BINDING - Python binding is requested
dnl  WITH_PYTHON_BINDING   - Either has swig or python wraper files
dnl ========================================================================
# Checks for Python (HAVE_PYTHON)
AM_PATH_PYTHON([2.2],
               [],
               [AC_MSG_WARN([Python not found. Python is required to build python binding. Python can be obtained from http://www.python.org])])
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")

# Checks for SWIG (HAVE_SWIG)
AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")

# Python binding enabled?
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding (default=yes)]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
AM_CONDITIONAL(ENABLE_PYTHON_BINDING, test "x$ac_enable_python_binding" = "xyes")

# Checks for gammalib_wrap.cpp and gammalib.py
AC_CHECK_FILES([pyext/gammalib_wrap.cpp pyext/gammalib.py],,)
has_python_wrap="no"
if test "x$ac_cv_file_pyext_gammalib_wrap_cpp" = "xyes"; then
  if test "x$ac_cv_file_pyext_gammalib_py" = "xyes"; then
    has_python_wrap="yes"
  fi
fi
AM_CONDITIONAL(WITH_PYTHON_BINDING, test "x$has_python_wrap" = "xyes" -o "$SWIG")

# Inform about missing components and set host dependent variables
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    if test "x$has_python_wrap" != "xyes"; then
      if test ! "$SWIG"; then
        AC_MSG_WARN([Neither pyext/gammalib_wrap.cpp and pyext/gammalib.py nor swig has not been found])
        AC_MSG_WARN([Python binding for gammalib cannot be built])
        AC_MSG_WARN([Either use a gammalib distribution with gammalib_wrap.cpp and gammalib.py shipped or install swig on your system])
        AC_MSG_WARN([swig can be obtained from http://www.swig.org/])
      fi
    fi

    # Set host specific build prefix (needed for Solaris pycc bug)
    PYTHON_BUILD_PREFIX=
    case $host in
      sparc*-sun-solaris*)  PYTHON_BUILD_PREFIX="CC=CC" ;;
      *)                    PYTHON_BUILD_PREFIX="" ;;
    esac                    
    AC_SUBST(PYTHON_BUILD_PREFIX)
    
    # Set host specific extra libraries. Add all host specific libraries here that
    # need to be linked in the python module
    PYTHON_EXTRA_LIBS=
    case $host in
      sparc*-sun-solaris*)  PYTHON_EXTRA_LIBS="Cstd" ;;
      *)                    PYTHON_EXTRA_LIBS="" ;;
    esac                    
    AC_SUBST(PYTHON_EXTRA_LIBS)

  else
    AC_MSG_WARN([Python binding for gammalib cannot be built as Python seems not to be installed on your system])
  fi
else
  AC_MSG_NOTICE([Python binding for gammalib will not be built])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi


dnl ==========================
dnl Checks for cfitsio library
dnl ==========================
# Set host specific cfitsio dependencies
cfitsio_deps=
case $host in
  sparc*-sun-solaris*)  cfitsio_deps="-lm -lsocket" ;;
  *)                    cfitsio_deps="-lm" ;;
esac                    

# Search for cfitsio shared library in a some standard directories
for i in "$gammalib_prefix" /usr /usr/local /opt/local ; do
  if test -f $i/lib/libcfitsio.so -o -f $i/lib/libcfitsio.dylib -a -f $i/include/fitsio.h; then
    LIBDIR_CFITSIO=$i/lib
    INCDIR_CFITSIO=$i/include
    break
  fi
  if test -f $i/lib64/libcfitsio.so -o -f $i/lib64/libcfitsio.dylib -a -f $i/include/fitsio.h; then
    LIBDIR_CFITSIO=$i/lib64
    INCDIR_CFITSIO=$i/include
    break
  fi
done
if test -n "$LIBDIR_CFITSIO" -a -n "$INCDIR_CFITSIO"; then
  has_cfitsio="yes"
  AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])
  LDFLAGS="$LDFLAGS -L$LIBDIR_CFITSIO"
  CPPFLAGS="$CPPFLAGS -I$INCDIR_CFITSIO"
  LIBS="${LIBS} -lcfitsio ${cfitsio_deps}"
  AC_SUBST(LIBDIR_CFITSIO)
  AC_SUBST(INCDIR_CFITSIO)
else
    has_cfitsio="no"
fi
AM_CONDITIONAL(WITH_CFITSIO, test "x$has_cfitsio" = "xyes")

# If cfitsio has not yet been found then search now using the AC_CHECK_LIB
# macro. This allows also search in paths that are eventually specified in
# LDFLAGS
if test "x$has_cfitsio" != "xyes"; then
  AC_CHECK_LIB([cfitsio], [ffpss],
               [LIBS="${LIBS} -lcfitsio ${cfitsio_deps}"
                has_cfitsio="yes"
                AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])],
               [], [${cfitsio_deps}])
  AM_CONDITIONAL(WITH_CFITSIO, test "x$ac_cv_lib_cfitsio_ffpss" = "xyes")
fi

## If cfitsio is available then get library path for unit testing
#if test "x$ac_cv_lib_cfitsio_ffpss" = "xyes"; then
#  for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
#    test -f $i/libcfitsio.a -o -f $i/libcfitsio.so -o -f $i/libcfitsio.dylib && LIBDIR_CFITSIO=$i
#  done
#  AC_SUBST(LIBDIR_CFITSIO)
#fi
#
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

# Search for readline shared library
for i in "$gammalib_prefix" /usr /usr/local /opt/local ; do
  if test -f $i/lib/libreadline.so -o -f $i/lib/libreadline.dylib -a -f $i/include/readline/readline.h; then
    LIBDIR_READLINE=$i/lib
    INCDIR_READLINE=$i/include
    break
  fi
  if test -f $i/lib64/libreadline.so -o -f $i/lib64/libreadline.dylib -a -f $i/include/readline/readline.h; then
    LIBDIR_READLINE=$i/lib64
    INCDIR_READLINE=$i/include
    break
  fi
done
if test -n "$LIBDIR_READLINE" -a -n "$INCDIR_READLINE"; then
  has_readline="yes"
  AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])
  if test "$LIBDIR_READLINE" != "$LIBDIR_CFITSIO"; then
    LDFLAGS="$LDFLAGS -L$LIBDIR_READLINE"
  fi
  if test "$INCDIR_READLINE" != "$INCDIR_CFITSIO"; then
    CPPFLAGS="$CPPFLAGS -I$INCDIR_READLINE"
  fi
  LIBS="${LIBS} -lreadline ${readline_deps}"
  AC_SUBST(LIBDIR_READLINE)
  AC_SUBST(INCDIR_READLINE)
else
    has_readline="no"
fi
AM_CONDITIONAL(WITH_READLINE, test "x$has_readline" = "xyes")

# If readline has not yet been found then search now using the AC_CHECK_LIB
# macro. This allows also search in paths that are eventually specified in
# LDFLAGS
if test "x$has_readline" != "xyes"; then
  AC_CHECK_LIB([readline], [main],
               [LIBS="${LIBS} -lreadline ${readline_deps}"
                has_readline="yes"
                AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])],
               [], [${readline_deps}])
  AM_CONDITIONAL(WITH_READLINE, test "x$ac_cv_lib_readline_main" = "xyes")
fi

# Check if readline is available
#AC_CHECK_LIB([readline], [main],
#             [LIBS="${LIBS} -lreadline ${readline_deps}"
#              AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])],
#             [], [${readline_deps}])
#AM_CONDITIONAL(WITH_READLINE, test "x$ac_cv_lib_readline_main" = "xyes")
#
# If readline is available then get library path for unit testing
#if test "x$ac_cv_lib_readline_main" = "xyes"; then
#  for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
#    test -f $i/libreadline.a -o -f $i/libreadline.so -o -f $i/libreadline.dylib && LIBDIR_READLINE=$i
#  done
#  AC_SUBST(LIBDIR_READLINE)
#fi
#
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


dnl =======
dnl Doxygen
dnl =======
DX_HTML_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_INIT_DOXYGEN([GammaLib], [doc/Doxyfile], [doc/doxygen])


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
                 src/gammalib-setup
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
if test "x$has_cfitsio" = "xyes"; then
  echo "  * FITS I/O support             (yes)   $LIBDIR_CFITSIO $INCDIR_CFITSIO"
else
  echo "  - FITS I/O support             (no)    no cfitsio library found"
fi
if test "x$has_readline" = "xyes"; then
  echo "  * Readline support             (yes)   $LIBDIR_READLINE $INCDIR_READLINE"
else
  echo "  - Readline support             (no)    no readline library found"
fi
if test ! "$SWIG"; then
  echo "  - swig                         (no)"
else
  echo "  * swig                         (yes)"
fi
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    if test "x$has_python_wrap" = "xyes"; then
      if test ! "$SWIG"; then
        echo "  * Python binding               (yes)"
      else
        echo "  * Python binding               (yes)   use swig for updates"
      fi
    else
      if test ! "$SWIG"; then
        echo "  - Python binding               (no)    swig not installed"
      else
        echo "  * Python binding               (yes)   use swig for building"
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
if test "$DX_DOXYGEN" != ""; then
  echo "  * Doxygen                      (yes)   $DX_DOXYGEN"
else
  echo "  - Doxygen                      (no)"
fi
echo
echo "Now type 'make'"
echo
