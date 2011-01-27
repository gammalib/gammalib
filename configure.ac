dnl =============================================================
dnl Process this file with autoconf to produce a configure script
dnl =============================================================
AC_INIT([GammaLib], [0.4], [knodlseder@cesr.fr])
AC_CONFIG_SRCDIR(src/Makefile.am)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu no-dependencies])
AC_CONFIG_HEADERS([config.h])
dnl AM_MAINTAINER_MODE


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
dnl AC_TYPE_SIZE_T
dnl AC_C_CONST
dnl AC_C_VOLATILE
dnl AC_C_INLINE
dnl gl_VISIBILITY


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


dnl ======================================================================
dnl Checks for python binding. Python binding requires Python and either a
dnl gammalib_wrap.cpp file or SWIG to create this file.
dnl ======================================================================
# Checks for SWIG (HAVE_SWIG)
AC_PATH_PROG([SWIG], [swig])
AM_CONDITIONAL(HAVE_SWIG, test "$SWIG")

# Checks for Python (HAVE_PYTHON)
AM_PATH_PYTHON([2.2],
               [],
               [AC_MSG_WARN([Python not found. Python is required to build python binding. Python can be obtained from http://www.python.org])])
AM_CONDITIONAL(HAVE_PYTHON, test "$PYTHON")

# Checks for gammalib_wrap.cpp and gammalib.py
AC_CHECK_FILES([pyext/gammalib_wrap.cpp pyext/gammalib.py],
               [ac_has_python_wrap="yes"],
               [ac_has_python_wrap="no"])
if test "x$ac_has_python_wrap" = "xyes"; then
  AC_MSG_NOTICE([Found gammalib_wrap.cpp and gammalib.py])
else
  AC_MSG_NOTICE([gammalib_wrap.cpp and/or gammalib.py not found])
fi
AM_CONDITIONAL(HAS_PYTHON_WRAP, test "x$ac_has_python_wrap" = "xyes")

# Python bindings
AC_ARG_ENABLE([python-binding],
              AS_HELP_STRING([--enable-python-binding],
                             [build python binding (default=yes)]),
              [ac_enable_python_binding="$enableval"],
              [ac_enable_python_binding="yes"])
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "$PYTHON"; then
    if test "x$ac_has_python_wrap" = "xyes"; then
      AC_MSG_NOTICE([Build Python binding using shipped gammalib_wrap.cpp and gammalib.py files])
    else
      if test ! "$SWIG"; then
        AC_MSG_WARN([Neither gammalib_wrap.cpp and/or gammalib.py files nor SWIG has not been found.])
        AC_MSG_WARN([Python binding for gammalib cannot be built.])
        AC_MSG_WARN([Either use a gammalib distribution with gammalib_wrap.cpp and gammalib.py shipped or install SWIG on your system.])
      fi
    fi
  else
    AC_MSG_WARN([Python binding for gammalib cannot be built as Python seems not to be installed on your system.])
  fi
else
  AC_MSG_NOTICE([Python binding for gammalib will not be built.])
  AC_MSG_NOTICE([Enable Python binding module building with --enable-python-binding])
fi
AM_CONDITIONAL(ENABLE_PYTHON_BINDING, test "x$ac_enable_python_binding" = "xyes")


dnl ==========================
dnl Checks for cfitsio library
dnl ==========================
#AC_SEARCH_LIBS(ffpss, cfitsio,
#               AC_DEFINE([G_CFITSIO], [1], [Define if cfitsio is available]), ,[-lm])
# Check if we want to use cfitsio
AC_ARG_WITH([cfitsio],
            [AS_HELP_STRING([--with-cfitsio],
                            [Use CFITSIO library (needed to support FITS files)])],
            [],
            [with_cfitsio=yes])
AC_ARG_WITH([cfitsio-lib-dir],
            AC_HELP_STRING([--with-cfitsio-lib-dir],
                           [CFITSIO library directory [[/usr/(local/)lib(64)]]]),
            [cfitsio_lib_dir=$withval])

# If cfitsio is requested then add cfitsio library path to LDFLAGS
if test "x$with_cfitsio" = "xyes"; then
  if test "x$cfitsio_lib_dir" != "x"; then
    test -f $cfitsio_lib_dir/libcfitsio.a -o -f $cfitsio_lib_dir/libcfitsio.so -o -f $cfitsio_lib_dir/libcfitsio.dylib && LIBDIR_CFITSIO=$cfitsio_lib_dir
  else
    for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
      test -f $i/libcfitsio.a -o -f $i/libcfitsio.so -o -f $i/libcfitsio.dylib && LIBDIR_CFITSIO=$i
    done
    for i in /usr/local/include /opt/local/include /usr/include ; do
      test -f $i/fitsio.h && INCDIR_CFITSIO=$i
    done
  fi
  if test -z "$LIBDIR_CFITSIO"; then
    if test "x$cfitsio_lib_dir" != "x";then
      AC_MSG_ERROR(Cannot find CFITSIO library under $cfitsio_lib_dir)
    else
      AC_MSG_ERROR(Cannot find CFITSIO library. Use --with-cfitsio-lib-dir= to specify non-default path or --without-cfitsio to compile with CFITSIO library.)
    fi
  fi
  LDFLAGS="$LDFLAGS -L$LIBDIR_CFITSIO"
  CPPFLAGS="$CPPFLAGS -I$LIBDIR_CFITSIO/../include"
  AC_SUBST(LIBDIR_CFITSIO)
fi

# Check for cfitsio library
if test "x$with_cfitsio" = "xyes"; then
  AC_CHECK_LIB([cfitsio], [ffpss],
               [LIBS="${LIBS} -lcfitsio -lm"
                AC_DEFINE([HAVE_LIBCFITSIO], [1], [Define if cfitsio library is available])],
               [AC_MSG_FAILURE([--with-cfitsio was given, but test for cfitsio failed])],
               [-lm])
fi
AM_CONDITIONAL(WITH_CFITSIO, test "x$with_cfitsio" = "xyes")


dnl ===========================
dnl Checks for readline library
dnl ===========================
# Checks for readline
AC_ARG_WITH([readline],
            [AS_HELP_STRING([--with-readline],
            [support fancy command line editing @<:@default=check@:>@])],
            [],
            [with_readline=yes])
AC_ARG_WITH([readline-lib-dir],
            AC_HELP_STRING([--with-readline-lib-dir],
                           [readline library directory [[/usr/(local/)lib(64)]]]),
            [readline_lib_dir=$withval])

# If readline is requested then add readline library path to LDFLAGS
if test "x$with_readline" = "xyes"; then
  if test "x$readline_lib_dir" != "x"; then
    test -f $readline_lib_dir/libreadline.a -o -f $readline_lib_dir/libreadline.so -o -f $readline_lib_dir/libreadline.dylib && LIBDIR_READLINE=$readline_lib_dir
  else
    for i in /usr/local/lib /usr/local/lib64 /opt/local/lib /usr/lib /usr/lib64 ; do
      test -f $i/libreadline.a -o -f $i/libreadline.so -o -f $i/libreadline.dylib && LIBDIR_READLINE=$i
    done
  fi
  if test -z "$LIBDIR_READLINE"; then
    if test "x$readline_lib_dir" != "x";then
      AC_MSG_ERROR(Cannot find readline library under $readline_lib_dir)
    else
      AC_MSG_ERROR(Cannot find readline library. Use --with-readline-lib-dir= to specify non-default path.)
    fi
  fi
  LDFLAGS="$LDFLAGS -L$LIBDIR_READLINE"
  AC_SUBST(LIBDIR_READLINE)
fi

# Check for readline library
if test "x$with_readline" = "xyes"; then
  AC_CHECK_LIB([readline], [main],
               [LIBS="${LIBS} -lreadline -lncurses"
                AC_DEFINE([HAVE_LIBREADLINE], [1], [Define if readline library is available])],
               [AC_MSG_FAILURE([--with-readline was given, but test for readline failed])],
               [-lncurses])
fi
AM_CONDITIONAL(WITH_READLINE, test "x$with_readline" = "xyes")

# Older readline.h's declared readline() without any arguments and outside of 'extern "C"', which freaks out C++
if test "$with_readline" != no; then
  AC_LANG_PUSH(C++)
  AC_MSG_CHECKING([if readline.h needs extern "C"])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stdio.h>
#include <readline/readline.h>
                                     ]],[[
  char*x=readline("test"); 
  return (x?0:1);
                                     ]])],
                    [AC_MSG_RESULT(no)], 
                    [AC_DEFINE(READLINE_H_NEEDS_EXTERN_C, 1, [Define to 1 if readline/readline.h needs extern "C"]) 
                     AC_MSG_RESULT(yes)])
  AC_LANG_POP(C++)
fi

# Slightly newer but still broken readline.h's use extern "C" but don't declare arguments as const
if test "$with_readline" != no; then
  AC_LANG_PUSH(C++)
  AC_MSG_CHECKING([if readline.h lacks 'const' in its declarations])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stdio.h>
#if READLINE_H_NEEDS_EXTERN_C
extern "C" {
#endif
#include <readline/readline.h>
#if READLINE_H_NEEDS_EXTERN_C
} // terminate extern "C"
#endif
                                     ]],[[
  const char* prompt = "abc";
  char*x = readline(prompt);
                                     ]])],
                     [AC_MSG_RESULT(no)],
                     [AC_DEFINE(READLINE_H_USES_NO_CONST, 1, [Define to 1 if readline/readline.h doesn't properly declare arguments to be const])
                      AC_MSG_RESULT(yes)])
  AC_LANG_POP(C++)
fi

# Some other broken readline.h's use extern "C" and declare arguments for most everything, except for the callback functions
if test "$with_readline" != no; then
  AC_LANG_PUSH(C++)
  AC_MSG_CHECKING([if readline.h lacks types in its declaration of callbacks])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <stdio.h>
#if READLINE_H_NEEDS_EXTERN_C
extern "C" {
#endif
#include <readline/readline.h>
#if READLINE_H_NEEDS_EXTERN_C
} // terminate extern "C"
#endif

#ifdef READLINE_H_USES_NO_CONST
char* dummy_completion(char*, int) { return 0; }
#else
char* dummy_completion(const char*, int) { return 0; }
#endif
                                     ]],[[
  rl_completion_entry_function = dummy_completion;
                                     ]])],
                     [AC_MSG_RESULT(no)],
                     [AC_DEFINE(READLINE_H_LACKS_TYPES_FOR_CALLBACKS, 1, [Define to 1 if readline/readline.h lacks proper declarations for callbacks])
                      AC_MSG_RESULT(yes)])
  AC_LANG_POP(C++)
fi


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
