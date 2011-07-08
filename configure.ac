# Configuration for GammaLib - a versatile toolbox for high-level analysis
#                              of astronomical gamma-ray data
#
# Copyright (C) 2006-2011 Jurgen Knodlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Process this file with autoconf to produce a configure script

AC_INIT([GammaLib], [0.4.9], [knodlseder@cesr.fr], [gammalib])
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
  gammalib_prefix=${prefix}
fi
AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${gammalib_prefix}", [Location where library is installed])


dnl ===========================
dnl Things required by automake
dnl ===========================
# AC_ARG_PROGRAM
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
# AC_TYPE_LONG_LONG_INT
AC_C_CONST
AC_C_INLINE


dnl ============================
dnl Checks for library functions
dnl ============================
# AC_CHECK_LIB(m, sincos, [AC_DEFINE([HAVE_SINCOS],[1],[Define to 1 if your system has `sincos'.])])


dnl ===================
dnl Checks for programs
dnl ===================
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL


dnl =========================================
dnl Checks for debugging or profiling options
dnl =========================================

# Check for debugging option
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE([debug],
              AS_HELP_STRING([--enable-debug],
                             [turn on debugging [default=no]]),
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

# Check for profiling option
AC_MSG_CHECKING(whether to enable profiling)
AC_ARG_ENABLE([profiling],
              AS_HELP_STRING([--enable-profiling],
                             [turn on profiling [default=no]]),
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
dnl NaN/Inf checking option (G_NAN_CHECK)
dnl =====================================
AC_MSG_CHECKING(whether to enable NaN/Inf checking)
AC_ARG_ENABLE([nan_check],
              AS_HELP_STRING([--enable-nan-check],
                             [performs NaN/Inf checking [default=yes]]),
             [NAN_CHECK="$enableval"],
             [NAN_CHECK="yes"])
if test "x$NAN_CHECK" = "xyes"; then
  AC_DEFINE([G_NAN_CHECK], [1], [Define if NaN/Inf checking should be performed])
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
                             [performs index range checking [default=yes]]),
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
                             [optimizes code for small memory [default=yes]]),
             [SMALL_MEMORY="$enableval"],
             [SMALL_MEMORY="yes"])
if test "x$SMALL_MEMORY" = "xyes"; then
  AC_DEFINE([G_SMALL_MEMORY], [1], [Define if small memory optimization is selected])
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi


dnl =========================================================================
dnl Checks if Python bindings are to be build. Python binding requires Python
dnl and either the files gammalib_wrap.cpp and gammalib.py or swig to create
dnl these files. The following conditionals are set:
dnl  WITH_PYTHON - Python bindings will be built
dnl  HAVE_SWIG    - swig is available
dnl =========================================================================

# Initialise flags
enable_python="no"
has_python="no"
has_python_h="no"
has_swig="no"
has_wrappers="no"

# Python binding enabled?
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding [default=yes]]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
if test "x$ac_enable_python_binding" = "xno"; then
  AC_MSG_NOTICE([Python binding for gammalib will not be built])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi

# Do we have Python and Python.h?
AM_PATH_PYTHON([2.2],
               [has_python="yes"],
               [AC_MSG_WARN([Python >= 2.2 not found. Python >= 2.2 is required to build python binding. Python can be obtained from http://www.python.org])])
if test "x$has_python" = "xyes"; then

  # Determine include path for python
  AC_PATH_PROG(PY_INCLUDE, python)  
  PY_INCLUDE=`echo $PY_INCLUDE | sed -e "s/bin/include/"`
  PY_INCLUDE="$PY_INCLUDE$PYTHON_VERSION"
  
  # Check if we have Python.h
  AC_CHECK_HEADERS($PY_INCLUDE/Python.h, [has_python_h="yes"], [has_python_h="no"])
  
  # Mac OS X kluge: if test before was unsuccessful then check for headers in /usr/include
  if test "x$has_python_h" = "xno"; then
    AC_CHECK_HEADERS("/usr/include/python$PYTHON_VERSION/Python.h", [has_python_h="yes"], [has_python_h="no"])
  fi

  # Signal any problem
  if test "x$has_python_h" = "xno"; then
    AC_MSG_WARN([No Python header file "$PY_INCLUDE/Python.h" found. Please install appropriate Python development package.])
  fi
fi

# Do we have swig?
AC_PATH_PROG([SWIG], [swig])
if test "$SWIG"; then
  has_swig="yes"
fi
AM_CONDITIONAL(HAVE_SWIG, test "x$has_swig" = "xyes")
 
# Do we have wrappers (gammalib_wrap.cpp and gammalib.py)?
AC_CHECK_FILES([pyext/gammalib_wrap.cpp pyext/gammalib.py],,)
if test "x$ac_cv_file_pyext_gammalib_wrap_cpp" = "xyes"; then
  if test "x$ac_cv_file_pyext_gammalib_py" = "xyes"; then
    has_wrappers="yes"
  fi
fi

# Signal if Python is enabled
if test "x$ac_enable_python_binding" = "xyes" -a "x$has_python" = "xyes" -a "x$has_python_h" = "xyes"; then
  if test "x$has_swig" = "xyes" -o "x$has_wrappers" = "xyes"; then
    enable_python="yes"
  fi
fi
AM_CONDITIONAL(WITH_PYTHON, test "x$enable_python" = "xyes")

# Set Python environment
if test "x$enable_python" = "xyes"; then

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
  if test "x$has_python_wrap" = "xno" -a "x$has_swig" = "xno"; then
    AC_MSG_WARN([Neither pyext/gammalib_wrap.cpp and pyext/gammalib.py nor swig has not been found])
    AC_MSG_WARN([Python binding for gammalib cannot be built])
    AC_MSG_WARN([Either use a gammalib distribution with gammalib_wrap.cpp and gammalib.py shipped or install swig on your system])
    AC_MSG_WARN([swig can be obtained from http://www.swig.org/])
  fi
fi


dnl ==========================
dnl Checks for cfitsio library
dnl ==========================

# Initialise cfitsio flag
has_cfitsio="no"

# Check if we want to use cfitsio
AC_ARG_WITH([cfitsio],
            [AS_HELP_STRING([--with-cfitsio],
                            [Use CFITSIO library [default=yes]])],
            [],
            [with_cfitsio=yes])

# If we want to use cfitsio, then search the library now
if test "x$with_cfitsio" = "xyes"; then

  # Set host specific cfitsio dependencies
  cfitsio_deps=
  case $host in
    sparc*-sun-solaris*)  cfitsio_deps="-lm -lsocket" ;;
    *)                    cfitsio_deps="-lm" ;;
  esac
  
  # First search for cfitsio in the standard paths using the AC_CHECK_LIB
  # macro. This allows also search in paths that are eventually specified by LDFAGS.
  AC_CHECK_LIB([cfitsio], [ffpss],
               [has_cfitsio="yes"
                AC_MSG_NOTICE([cfitsio found])],
               [], [${cfitsio_deps}])

  # If cfitsio has not been found then search for cfitsio shared library in a some
  # default directories
  if test "x$has_cfitsio" = "xno"; then
    for i in "$gammalib_prefix" /opt/local /usr/local /usr ; do
      if test -f $i/lib64/libcfitsio.so -o -f $i/lib64/libcfitsio.dylib -a -f $i/include/fitsio.h; then
        LIBDIR_CFITSIO=$i/lib64
        INCDIR_CFITSIO=$i/include
        break
      fi
      if test -f $i/lib/libcfitsio.so -o -f $i/lib/libcfitsio.dylib -a -f $i/include/fitsio.h; then
        LIBDIR_CFITSIO=$i/lib
        INCDIR_CFITSIO=$i/include
        break
      fi
    done
    if test -n "$LIBDIR_CFITSIO" -a -n "$INCDIR_CFITSIO"; then
      has_cfitsio="yes"
      LDFLAGS="$LDFLAGS -L$LIBDIR_CFITSIO"
      CPPFLAGS="$CPPFLAGS -I$INCDIR_CFITSIO"
      AC_MSG_NOTICE([cfitsio found in ${LIBDIR_CFITSIO}])
    fi
  fi

  # Set cfitsio information
  if test "x$has_cfitsio" = "xyes"; then
    LIBS_CFITSIO="-lcfitsio ${cfitsio_deps}"
    LIBS="${LIBS} ${LIBS_CFITSIO}"
  fi

fi
if test "x$has_cfitsio" = "xyes"; then
  AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])
  AC_SUBST(CFITSIO)
  AC_SUBST(LIBS_CFITSIO)
  AC_SUBST(LIBDIR_CFITSIO)
  AC_SUBST(INCDIR_CFITSIO)
fi
AM_CONDITIONAL(WITH_CFITSIO, test "x$has_cfitsio" = "xyes")


dnl ===========================
dnl Checks for readline library
dnl ===========================

# Initialise readline flag
has_readline="no"

# Check if we want to use readline
AC_ARG_WITH([readline],
            [AS_HELP_STRING([--with-readline],
                            [Use readline library [default=yes]])],
            [],
            [with_readline=yes])

# If we want to use readline, then search now the library
if test "x$with_readline" = "xyes"; then

  # Set host specific readline dependencies
  readline_deps=
  case $host in
    sparc*-sun-solaris*)  readline_deps="-lncurses" ;;
    *)                    readline_deps="-lncurses" ;;
  esac
  
  # First search for readline using the AC_CHECK_LIB macro. This performs search
  # in all standard locations and the paths that are eventually specified by LDFLAGS
  AC_CHECK_LIB([readline], [main],
               [has_readline="yes"
                AC_MSG_NOTICE([readline found])],
               [], [${readline_deps}])

  # If readline has not been found then search for readline shared library in a couple
  # of predefined locations
  if test "x$has_readline" = "xno"; then
    for i in "$gammalib_prefix" /opt/local /usr/local /usr ; do
      if test -f $i/lib64/libreadline.so -o -f $i/lib64/libreadline.dylib -a -f $i/include/readline/readline.h; then
        LIBDIR_READLINE=$i/lib64
        INCDIR_READLINE=$i/include
        break
      fi
      if test -f $i/lib/libreadline.so -o -f $i/lib/libreadline.dylib -a -f $i/include/readline/readline.h; then
        LIBDIR_READLINE=$i/lib
        INCDIR_READLINE=$i/include
        break
      fi
    done
    if test -n "$LIBDIR_READLINE" -a -n "$INCDIR_READLINE"; then
      has_readline="yes"
      if test "$LIBDIR_READLINE" != "$LIBDIR_CFITSIO"; then
        LDFLAGS="$LDFLAGS -L$LIBDIR_READLINE"
      fi
      if test "$INCDIR_READLINE" != "$INCDIR_CFITSIO"; then
        CPPFLAGS="$CPPFLAGS -I$INCDIR_READLINE"
      fi
      AC_SUBST(LIBDIR_READLINE)
      AC_SUBST(INCDIR_READLINE)
      AC_MSG_NOTICE([readline found in ${LIBDIR_READLINE}])
    fi
  fi

  # Set readline information
  if test "x$has_readline" = "xyes"; then
    LIBS_READLINE="-lreadline ${readline_deps}"
    AC_SUBST(LIBS_READLINE)
    LIBS="${LIBS} ${LIBS_READLINE}"
  fi

fi
if test "x$has_readline" = "xyes"; then
  AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])
fi
AM_CONDITIONAL(WITH_READLINE, test "x$has_readline" = "xyes")


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
if test "x$has_python" = "xyes"; then
  echo "  * Python                       (yes)"
else
  echo "  - Python                       (no)"
fi
if test "x$has_python_h" = "xyes"; then
  echo "  * Python.h                     (yes)"
else
  echo "  - Python.h                     (no)"
fi
if test "x$has_swig" = "xyes"; then
  echo "  * swig                         (yes)"
else
  echo "  - swig                         (no)"
fi
if test "x$has_wrappers" = "xyes"; then
  echo "  * Python wrappers              (yes)"
else
  echo "  - Python wrappers              (no)"
fi

# Dump Python binding information
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "x$has_python" = "xyes"; then
    if test "x$has_python_h" = "xyes"; then
      if test "x$has_wrappers" = "xyes"; then
        if test "x$has_swig" = "xyes"; then
          echo "  * Make Python binding          (yes)   use swig for updates"
        else
          echo "  * Make Python binding          (yes)   use wrappers"
        fi
      else
        if test "x$has_swig" = "xyes"; then
          echo "  * Make Python binding          (yes)   use swig for building"
        else
          echo "  - Make Python binding          (no)    swig required to build wrappers"
        fi
      fi
    else
      echo "  - Make Python binding          (no)    Python.h not installed"
    fi
  else
    echo "  - Make Python binding          (no)    Python not installed"
  fi
else
  echo "  - Make Python binding          (no)"
fi

# Dump instrument interface information
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

# Dump Doxygen information
if test "$DX_DOXYGEN" != ""; then
  echo "  * Doxygen                      (yes)   $DX_DOXYGEN"
else
  echo "  - Doxygen                      (no)"
fi

# Dump compile options
if test "x$NAN_CHECK" = "xyes"; then
  echo "  * Perform in NaN/Inf checks    (yes)   (default)"
else
  echo "  - Perform in NaN/Inf checks    (no)"
fi
if test "x$RANGE_CHECK" = "xyes"; then
  echo "  * Perform range checking       (yes)   (default)"
else
  echo "  - Perform range checking       (no)"
fi
if test "x$SMALL_MEMORY" = "xyes"; then
  echo "  * Optimize memory usage        (yes)   (default)"
else
  echo "  - Optimize memory usage        (no)"
fi
if test "x$enable_debug" = "xyes"; then
  echo "  * Compile in debug code        (yes)"
else
  echo "  - Compile in debug code        (no)    (default)"
fi
if test "x$enable_profiling" = "xyes"; then
  echo "  * Enable code for profiling    (yes)"
else
  echo "  - Enable code for profiling    (no)    (default)"
fi

# Dump next step information
echo
echo "Now type 'make'"
echo
