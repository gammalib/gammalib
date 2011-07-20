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

AC_INIT([GammaLib], [0.4.10], [knodlseder@cesr.fr], [gammalib])
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


dnl ==============================================================
dnl Checks for programs
dnl (this has to come before the check for header files, otherwise
dnl  we won't detect if the C++ compiler is missing)
dnl ==============================================================
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL


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

  # Do we have distutils?
  AC_MSG_CHECKING([for the distutils Python package])
  ac_distutils_result=`$PYTHON -c "import distutils" 2>&1`
  if test -z "$ac_distutils_result"; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
  fi

  # If we have distutils then get python path from distutils ...
  if test -z "$ac_distutils_result"; then
    AC_MSG_CHECKING([for Python include path])
    python_path=`$PYTHON -c "import distutils.sysconfig; \
           		             print distutils.sysconfig.get_python_inc();"`
    if test -n "${python_path}"; then
      python_path="-I$python_path"
      has_python_h="yes"
    fi
    PYTHON_CPPFLAGS=$python_path
    CPPFLAGS="${CPPFLAGS} ${python_path}"
	AC_MSG_RESULT([$python_path])
	AC_SUBST([PYTHON_CPPFLAGS])

  # ... otherwise check for includes relative to exectuable (obsolete?)
  else                             
                             
    # Determine include path for python
    AC_PATH_PROG(PY_INCLUDE, python)  
    PY_INCLUDE=`echo $PY_INCLUDE | sed -e "s/bin/include/"`
    PY_INCLUDE="$PY_INCLUDE$PYTHON_VERSION"
  
    # Check if we have Python.h
    AC_CHECK_HEADERS($PY_INCLUDE/Python.h, [has_python_h="yes"], [has_python_h="no"])
  
    # Signal any problem
    if test "x$has_python_h" = "xno"; then
      AC_MSG_WARN([No Python header file "$PY_INCLUDE/Python.h" found. Please install appropriate Python development package.])
    fi
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


dnl ===========================
dnl Checks for readline library
dnl ===========================

# Initialise readline and ncurses flags
has_readline="no"
has_ncurses="no"

# Check if we want to use readline
AC_ARG_WITH([readline],
            [AS_HELP_STRING([--with-readline],
                            [Use readline library [default=yes]])],
            [],
            [with_readline=yes])

# If we want to use readline, then search now the library
if test "x$with_readline" = "xyes"; then

  # First check for the presence of the ncurses library
  AC_FIND_LIB_HEADER([ncurses], [tgetstr], [], [has_ncurses="yes"])
  AC_FIND_LIB_HEADER([readline], [readline], [readline/readline.h], [has_readline="yes"])

fi
if test "x$has_readline" = "xyes"; then
  AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])
fi
AM_CONDITIONAL(WITH_READLINE, test "x$has_readline" = "xyes")


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

  # Is pkgconfig available? This sets PKG_CONFIG to the binary location
  PKG_PROG_PKG_CONFIG

  # Use pkgconfig if available
  if test "x$PKG_CONFIG" != "x"; then
    PKG_CHECK_MODULES([cfitsio], [cfitsio],
                      [has_cfitsio="yes"],
                      [AC_MSG_WARN([cfitsio pkgconfig information not found])])
    if test "x$has_cfitsio" = "xyes"; then
    
      # OpenSUSE kluge: cfitsio pkgconfig not correct
      if test -f "/usr/include/libcfitsio0/fitsio.h"; then
        cfitsio_CFLAGS="-I/usr/include/libcfitsio0"
      fi
      
      # Set environment variables
      CFLAGS="$CFLAGS $cfitsio_CFLAGS"
      CPPFLAGS="$CPPFLAGS $cfitsio_CFLAGS"
      LDFLAGS="$LDFLAGS $cfitsio_LIBS"
      CFITSIO="cfitsio"
    fi
  fi
    
  # If we still don't have cfitsio, then search the library and headers in
  # standard paths
  if test "x$has_cfitsio" = "xno"; then

    # Set host specific cfitsio dependencies
    cfitsio_deps=
    case $host in
      *-*-solaris*)  cfitsio_deps="-lm -lsocket" ;;
      *)             cfitsio_deps="-lm" ;;
    esac
  
    # Check for cfitsio library and headers. The cfitsio/fitsio.h is for
    # Scientific Linux as no pkgconfig is installed for this distribution
    AC_FIND_LIB_HEADER([cfitsio], [ffpss], [fitsio.h cfitsio/fitsio.h],
                       [has_cfitsio="yes"], [], [${cfitsio_deps}])

    # Add cfitsio dependencies
    if test "x$has_cfitsio" = "xyes"; then
      LIBS="${LIBS} ${cfitsio_deps}"
    fi
    
  fi
  
fi
if test "x$has_cfitsio" = "xyes"; then
  AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])
fi
AM_CONDITIONAL(WITH_CFITSIO, test "x$has_cfitsio" = "xyes")


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
  echo "  * FITS I/O support             (yes)"
else
  echo "  - FITS I/O support             (no)    no cfitsio library found"
fi
if test "x$has_readline" = "xyes"; then
  echo "  * Readline support             (yes)"
else
  if test "x$has_ncurses" = "xyes"; then
    echo "  - Readline support             (no)    no readline library found"
  else
    echo "  - Readline support             (no)    no ncurses library found"
  fi
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
