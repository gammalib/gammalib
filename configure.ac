dnl =============================================================
dnl Process this file with autoconf to produce a configure script
dnl =============================================================
AC_INIT([GammaLib], [0.4], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR(src/Makefile.am)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu no-dependencies])
AC_CONFIG_HEADERS([config.h])
dnl AM_MAINTAINER_MODE


# ==========
# Set prefix
# ==========
AC_PREFIX_DEFAULT([/usr/local/gamma])
if test "x$prefix" = xNONE; then
  AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], ["/usr/local/gamma"], 
		     [Location where library is installed])
else
  AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${prefix}", 
  		     [Location where library is installed])
fi


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


dnl ===========================
dnl things required by automake
dnl ===========================
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET


dnl ===================
dnl Checks for programs
dnl ===================
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Checks for SWIG
AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")

# Checks for Python
AM_PATH_PYTHON([2.2],
	       [],
	       [AC_MSG_WARN([Python not found. Python is required to build python binding. Python can be obtained from http://www.python.org])])
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")

AC_ARG_ENABLE([python-binding],
	      AS_HELP_STRING([--enable-python-binding],
			     [build python binding (default=yes)]),
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

# Checks for CTA interface
AC_ARG_WITH([cta],
            [AS_HELP_STRING([--with-cta],
                            [compile in CTA specific interface])],
            [],
            [with_cta=yes])
AM_CONDITIONAL(WITH_INST_CTA, test "x$with_cta" = "xyes")

# Checks for LAT interface
AC_ARG_WITH([lat],
            [AS_HELP_STRING([--with-lat],
                            [compile in LAT specific interface])],
            [],
            [with_lat=yes])
AM_CONDITIONAL(WITH_INST_LAT, test "x$with_lat" = "xyes")

# Checks for readline
#AC_ARG_WITH([readline],
#            [AS_HELP_STRING([--with-readline],
#            [support fancy command line editing @<:@default=check@:>@])],
#            [],
#            [with_readline=check])
#	  
#LIBREADLINE=
#AS_IF([test "x$with_readline" != xno],
#  [AC_CHECK_LIB([readline], [main],
#    [AC_SUBST([LIBREADLINE], ["-lreadline -lncurses"])
#     AC_DEFINE([HAVE_LIBREADLINE], [1],
#               [Define if you have libreadline])
#    ],
#    [if test "x$with_readline" != xcheck; then
#       AC_MSG_FAILURE(
#         [--with-readline was given, but test for readline failed])
#     fi
#    ], -lncurses)])


dnl ==========================
dnl Checks for cfitsio library
dnl ==========================
#AC_SEARCH_LIBS(ffpss, cfitsio,
#    AC_DEFINE([G_CFITSIO], [1], [Define if cfitsio is available]), ,[-lm])
# Check if we want to use cfitsio
CFITSIO="yes"
AC_ARG_WITH([cfitsio],
	        AC_HELP_STRING([--with-cfitsio],
		                   [Use CFITSIO]),
	        [CFITSIO=$withval])
AC_ARG_WITH([cfitsio-lib-dir],
	        AC_HELP_STRING([--with-cfitsio-lib-dir],
		                   [CFITSIO library directory [[/usr/(local/)lib[64]]]]),
	        [cfitsio_lib_dir=$withval])
if test "x$CFITSIO" = "xyes"; then
    AC_DEFINE([G_CFITSIO], [1], [Define if cfitsio is available])
fi

# If cfitsio is requested then add cfitsio library path to LDFLAGS
if test "x$CFITSIO" = "xyes"; then
  if test "x$cfitsio_lib_dir" != "x"; then
    test -f $cfitsio_lib_dir/libcfitsio.a -o -f $cfitsio_lib_dir/libcfitsio.so -o -f $cfitsio_lib_dir/libcfitsio.dylib && CFITSIO_LIBDIR=$cfitsio_lib_dir
  else
    for i in /usr/local/lib /usr/local/lib64 /usr/lib /usr/lib64 ; do
      test -f $i/libcfitsio.a -o -f $i/libcfitsio.so -o -f $i/libcfitsio.dylib && CFITSIO_LIBDIR=$i
    done
  fi
  if test -z "$CFITSIO_LIBDIR"; then
    if test "x$cfitsio_lib_dir" != "x";then
      AC_MSG_ERROR(Cannot find CFITSIO library under $cfitsio_lib_dir)
    else
      AC_MSG_ERROR(Cannot find CFITSIO library.  Use --with-cfitsio-lib-dir= to specify non-default path.)
    fi
  fi
  LDFLAGS="$LDFLAGS -L$CFITSIO_LIBDIR"
fi


dnl =======================
dnl Checks for header files
dnl =======================
AC_STDC_HEADERS


dnl =============================================================
dnl Checks for typedefs, structures, and compiler characteristics
dnl =============================================================
dnl AC_TYPE_SIZE_T
dnl AC_C_CONST
dnl AC_C_VOLATILE
dnl AC_C_INLINE
dnl gl_VISIBILITY


dnl ============================
dnl Checks for library functions
dnl ============================


dnl =========================
dnl Check for system services
dnl =========================
AC_CANONICAL_HOST


dnl =====================
dnl Range checking option
dnl =====================
RANGE_CHECK="yes"
AC_ARG_ENABLE(range_check,
             [  --enable-range-check    performs index range checking (default=yes)],
             [ RANGE_CHECK="$enableval" ])
if test "x$RANGE_CHECK" = "xyes"; then
        AC_DEFINE([G_RANGE_CHECK], [1], [Define if range checking should be performed])
fi


dnl ===================
dnl Small memory option
dnl ===================
SMALL_MEMORY="yes"
AC_ARG_ENABLE(small_memory,
             [  --enable-small-memory   optimizes code for small memory (default=yes)],
             [ SMALL_MEMORY="$enableval" ])
if test "x$SMALL_MEMORY" = "xyes"; then
        AC_DEFINE([G_SMALL_MEMORY], [1], [Define if small memory optimization is selected])
fi


dnl ======================================
dnl Set list of output files to be creates
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
                 inst/Makefile
                 include/Makefile
                 test/Makefile
                 pyext/Makefile
                 pyext/setup.py
                 gammalib.pc])
if test "x$with_cta" = "xyes"; then
    AC_CONFIG_FILES([inst/cta/Makefile])
fi
if test "x$with_lat" = "xyes"; then
    AC_CONFIG_FILES([inst/lat/Makefile])
fi
AC_OUTPUT
