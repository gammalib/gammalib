#############################################################################
# Configuration for GammaLib - a versatile toolbox for high-level analysis  #
#                              of astronomical gamma-ray data               #
# ------------------------------------------------------------------------- #
# Copyright (C) 2006-2011 Jurgen Knodlseder                                 #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
# ------------------------------------------------------------------------- #
# Process this file with autoconf to produce a configure script.            #
# ------------------------------------------------------------------------- #
# ToDo:                                                                     #
# - add possibility to build as Framework on Darwin                         #
# - strip "-Wstrict-prototypes" for C++ compiling                           #
#############################################################################

#############################################################################
# Initialisation                                                            #
#############################################################################
AC_INIT([gammalib], [0.5.0], [knodlseder@cesr.fr], [gammalib])
AC_CONFIG_SRCDIR([src/support/GTools.cpp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([gnu no-dependencies])
AC_CONFIG_HEADERS([config.h])


#############################################################################
# Library versioning                                                        #
# ------------------------------------------------------------------------- #
# See the libtool manual for an explanation.  Currently, library versioning #
# is not yet used for GammaLib as we're still in the development cycle.     #
# Users have to make sure that they link against the correct library.       #
#############################################################################
GAMMALIB_LT_VERSION="0:0:0"
AC_SUBST(GAMMALIB_LT_VERSION)


#############################################################################
# Set prefix                                                                #
# ------------------------------------------------------------------------- #
# By default, GammaLib gets installed into "/usr/local/gamma".  The default #
# can be overwritten using the --prefix option.                             #
#############################################################################
AC_PREFIX_DEFAULT([/usr/local/gamma])
if test "x$prefix" = xNONE; then
  gammalib_prefix="/usr/local/gamma"
else
  gammalib_prefix=${prefix}
fi
AC_DEFINE_UNQUOTED([PACKAGE_PREFIX], "${gammalib_prefix}", [Location where library is installed])


#############################################################################
# Check for programs that are required by automake                          #
#############################################################################
AC_PROG_MAKE_SET


#############################################################################
# Determine the build and host system type                                  #
#############################################################################
AC_CANONICAL_HOST


#############################################################################
# Determine the C and C++ compilers.  This has to come before the check for #
# header files, otherwise we won't detect if the C++ compiler is missing.   #
#############################################################################
AC_PROG_CC
AC_PROG_CXX


#############################################################################
# Darwin option: Determine architecture of universal library to be built    #
# ------------------------------------------------------------------------- #
# --with-univeral-archs=VALUE, where VALUE is one of "32-bit", "64-bit",    #
# "3-way", "intel" or "all".  Defaults to "32-bit".                         #
#############################################################################
UNIVERSAL_ARCHS="32-bit"
AC_SUBST(LIPO_32BIT_FLAGS)
AC_MSG_CHECKING(for --with-universal-archs)
AC_ARG_WITH(universal-archs,
    AS_HELP_STRING([--with-universal-archs=ARCH], [select architectures for universal build ("32-bit", "64-bit", "3-way", "intel" or "all")]),
[
	AC_MSG_RESULT($withval)
	UNIVERSAL_ARCHS="$withval"
],
[
 	AC_MSG_RESULT(32-bit)
])


#############################################################################
# Darwin option: Determine if a universal library should be built           #
# ------------------------------------------------------------------------- #
# --enable-universalsdk=PATH takes an optional argument that specifies      #
# which OSX SDK should be used to perform the build.                        #
# This defaults to:                                                         #
# "/Developer/SDKs/MacOSX.10.4u.sdk" for a 32-bit build (default)           #
# "/Developer/SDKs/MacOSX.10.5.sdk"  for any other build                    #
#############################################################################
AC_MSG_CHECKING([for --enable-universalsdk])
AC_ARG_ENABLE(universalsdk,
	AS_HELP_STRING([--enable-universalsdk@<:@=SDKDIR@:>@], [Build against Mac OS X 10.4u SDK (ppc/i386)]),
[
	case $enableval in
	yes)
        if [test "${UNIVERSAL_ARCHS}" = "32-bit";] then
            enableval=/Developer/SDKs/MacOSX10.4u.sdk
        else
            enableval=/Developer/SDKs/MacOSX10.5.sdk
        fi
		if [test ! -d "${enableval}";] then
			enableval=/
        fi
		;;
	esac
	case $enableval in
	no)
		UNIVERSALSDK=
		enable_universalsdk=
		;;
	*)
		UNIVERSALSDK=$enableval
		if [test ! -d "${UNIVERSALSDK}"]; then
			AC_MSG_ERROR([--enable-universalsdk specifies non-existing SDK: ${UNIVERSALSDK}])
		fi
		;;
	esac
	
],[
   	UNIVERSALSDK=
	enable_universalsdk=
])
if [test -n "${UNIVERSALSDK}"]; then
	AC_MSG_RESULT(${UNIVERSALSDK})
else
	AC_MSG_RESULT(no)
fi
AC_SUBST(UNIVERSALSDK)


#############################################################################
# Darwin option: Check whether --with-universal-archs was specified without #
#                --enable-universalsdk                                      #
#############################################################################
if [test "${withval}"]; then
    if [test "${enable_universalsdk}"]; then
		:
	else
		AC_MSG_ERROR([--with-universal-archs without --enable-universalsdk.])
	fi
fi


#############################################################################
# Record MACOSX_DEPLOYMENT_TARGET                                           #
# ------------------------------------------------------------------------- #
# Record the configure-time value of MACOSX_DEPLOYMENT_TARGET, it may       #
# influence the way we can build extensions, so distutils needs to check it #
#############################################################################
AC_SUBST(CONFIGURE_MACOSX_DEPLOYMENT_TARGET)
AC_SUBST(EXPORT_MACOSX_DEPLOYMENT_TARGET)
CONFIGURE_MACOSX_DEPLOYMENT_TARGET=
EXPORT_MACOSX_DEPLOYMENT_TARGET='#'


#############################################################################
# Set name for machine-dependent library files                              #
#############################################################################
AC_SUBST(MACHDEP)
AC_MSG_CHECKING(MACHDEP)
if test -z "$MACHDEP"; then
    ac_sys_system=`uname -s`
	if test "$ac_sys_system" = "AIX" \
	-o "$ac_sys_system" = "UnixWare" -o "$ac_sys_system" = "OpenUNIX"; then
		ac_sys_release=`uname -v`
	else
		ac_sys_release=`uname -r`
	fi
	ac_md_system=`echo $ac_sys_system |
			   tr -d '[/ ]' | tr '[[A-Z]]' '[[a-z]]'`
	ac_md_release=`echo $ac_sys_release |
			   tr -d '[/ ]' | sed 's/^[[A-Z]]\.//' | sed 's/\..*//'`

    # Set MACHDEP
    MACHDEP="$ac_md_system$ac_md_release"
    case $MACHDEP in
    cygwin*) MACHDEP="cygwin";;
    darwin*) MACHDEP="darwin";;
    atheos*) MACHDEP="atheos";;
    irix646) MACHDEP="irix6";;
    '')	MACHDEP="unknown";;
    esac
fi
AC_MSG_CHECKING(machine type as reported by uname -m)
ac_sys_machine=`uname -m`
AC_MSG_RESULT($ac_sys_machine)


#############################################################################
# Set -arch flags for universal builds on OSX and tweak BASECFLAGS based on #
# compiler and platform (note that BASECFLAGS are not yet used).            #
#############################################################################
# Initialise the -arch flags for universal builds on OSX
UNIVERSAL_ARCH_FLAGS=
AC_SUBST(UNIVERSAL_ARCH_FLAGS)

# Tweak BASECFLAGS based on compiler and platform
case $GCC in
yes)

    # if using gcc on alpha, use -mieee to get (near) full IEEE 754
    # support.  Without this, treatment of subnormals doesn't follow
    # the standard.
    case $ac_sys_machine in
        alpha*)
            BASECFLAGS="$BASECFLAGS -mieee"
            ;;
    esac

    case $ac_sys_system in
	SCO_SV*)
	    BASECFLAGS="$BASECFLAGS -m486 -DSCO5"
	    ;;
	# is there any other compiler on Darwin besides gcc?
	Darwin*)
	    # -Wno-long-double, -no-cpp-precomp, and -mno-fused-madd
	    # used to be here, but non-Apple gcc doesn't accept them.
	    if test "${CC}" = gcc; then
            AC_MSG_CHECKING(which compiler should be used)
		    case "${UNIVERSALSDK}" in
		    */MacOSX10.4u.sdk)
			    # Build using 10.4 SDK, force usage of gcc when the 
			    # compiler is gcc, otherwise the user will get very
			    # confusing error messages when building on OSX 10.6
			    CC=gcc-4.0
			    CXX=g++-4.0
			    CPP=cpp-4.0
                ac_ct_CC="$CC"
                ac_ct_CXX="$CXX"
                ac_cv_prog_ac_ct_CC="$ac_ct_CC"
                ac_cv_prog_ac_ct_CXX="$ac_ct_CXX"
			    ;;
		    esac
		    AC_MSG_RESULT($CXX)
        fi

	    # Calculate the right deployment target for this build.
	    cur_target=`sw_vers -productVersion | sed 's/\(10\.[[0-9]]*\).*/\1/'`
	    if test ${cur_target} '>' 10.2; then
		    cur_target=10.3
		    if test ${enable_universalsdk}; then
			    if test "${UNIVERSAL_ARCHS}" = "all"; then
				    # Ensure that the default platform for a 
				    # 4-way universal build is OSX 10.5, 
				    # that's the first OS release where 
				    # 4-way builds make sense.
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "3-way"; then
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "intel"; then
				    cur_target='10.5'
			    elif test "${UNIVERSAL_ARCHS}" = "64-bit"; then
				    cur_target='10.5'
			    elif test `/usr/bin/arch` = "i386"; then
				    # On Intel macs default to a deployment
				    # target of 10.4, that's the first OSX
				    # release with Intel support.
				    cur_target="10.4"
			    fi
		    else
			    if test `/usr/bin/arch` = "i386"; then
				    # On Intel macs default to a deployment
				    # target of 10.4, that's the first OSX
				    # release with Intel support.
				    cur_target="10.4"
			    fi
		    fi
	    fi
	    CONFIGURE_MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET-${cur_target}}
        
        # Make sure that deployment target is >= SDK
        case "${UNIVERSALSDK}" in
        */MacOSX10.4u.sdk)
            if [test ${cur_target} '>' 10.4;] then
	            AC_MSG_ERROR([Invalid MACOSX_DEPLOYMENT_TARGET=$CONFIGURE_MACOSX_DEPLOYMENT_TARGET for SDK 10.4. Specify --enable-universalsdk=/Developer/SDKs/MacOSX${cur_target}.sdk])
            fi
            ;;
        */MacOSX10.5.sdk)
            if [test ${cur_target} '>' 10.5;] then
	            AC_MSG_ERROR([Invalid MACOSX_DEPLOYMENT_TARGET=$CONFIGURE_MACOSX_DEPLOYMENT_TARGET for SDK 10.5. Specify --enable-universalsdk=/Developer/SDKs/MacOSX${cur_target}.sdk])
            fi
            ;;
        */MacOSX10.6.sdk)
            if [test ${cur_target} '>' 10.6;] then
	            AC_MSG_ERROR([Invalid MACOSX_DEPLOYMENT_TARGET=$CONFIGURE_MACOSX_DEPLOYMENT_TARGET for SDK 10.6. Specify --enable-universalsdk=/Developer/SDKs/MacOSX${cur_target}.sdk])
            fi
            ;;
        esac
	    
	    # Make sure that MACOSX_DEPLOYMENT_TARGET is set in the 
	    # environment with a value that is the same as what we'll use
	    # in the Makefile to ensure that we'll get the same compiler
	    # environment during configure and build time.
	    MACOSX_DEPLOYMENT_TARGET="$CONFIGURE_MACOSX_DEPLOYMENT_TARGET"
	    export MACOSX_DEPLOYMENT_TARGET
	    EXPORT_MACOSX_DEPLOYMENT_TARGET=''

        # Set compile, link and pre-processor flags dependent on architecture
	    if test "${enable_universalsdk}"; then
		    UNIVERSAL_ARCH_FLAGS=""
	        if test "$UNIVERSAL_ARCHS" = "32-bit" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch ppc -arch i386"
		        ARCH_RUN_32BIT=""
		        LIPO_32BIT_FLAGS=""
            elif test "$UNIVERSAL_ARCHS" = "64-bit" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch ppc64 -arch x86_64"
		        LIPO_32BIT_FLAGS=""
		        ARCH_RUN_32BIT="true"
            elif test "$UNIVERSAL_ARCHS" = "all" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch i386 -arch ppc -arch ppc64 -arch x86_64"
		        LIPO_32BIT_FLAGS="-extract ppc7400 -extract i386"
		        ARCH_RUN_32BIT="/usr/bin/arch -i386 -ppc"
            elif test "$UNIVERSAL_ARCHS" = "intel" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch i386 -arch x86_64"
		        LIPO_32BIT_FLAGS="-extract i386"
		        ARCH_RUN_32BIT="/usr/bin/arch -i386"
            elif test "$UNIVERSAL_ARCHS" = "3-way" ; then
		        UNIVERSAL_ARCH_FLAGS="-arch i386 -arch ppc -arch x86_64"
		        LIPO_32BIT_FLAGS="-extract ppc7400 -extract i386"
		        ARCH_RUN_32BIT="/usr/bin/arch -i386 -ppc"
		    else
	            AC_MSG_ERROR([proper usage is --with-universal-arch=32-bit|64-bit|all|intel|3-way])
		    fi
		    CXXFLAGS="${UNIVERSAL_ARCH_FLAGS} ${CXXFLAGS}"
		    if test "${UNIVERSALSDK}" != "/"; then
			    CPPFLAGS="-isysroot ${UNIVERSALSDK} ${CPPFLAGS}"
			    LDFLAGS="-isysroot ${UNIVERSALSDK} ${LDFLAGS}"
			    CXXFLAGS="-isysroot ${UNIVERSALSDK} ${CXXFLAGS}"
		    fi
	    fi
	    ;;
	OSF*)
	    BASECFLAGS="$BASECFLAGS -mieee"
	    ;;
    esac
    ;;
*)
    case $ac_sys_system in
    OpenUNIX*|UnixWare*)
	BASECFLAGS="$BASECFLAGS -K pentium,host,inline,loop_unroll,alloca "
	;;
    OSF*)
	BASECFLAGS="$BASECFLAGS -ieee -std"
    	;;
    SCO_SV*)
	BASECFLAGS="$BASECFLAGS -belf -Ki486 -DSCO5"
	;;
    esac
    ;;
esac


#############################################################################
# Checks for debugging and profiling options                                #
#############################################################################
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


#############################################################################
# Set linker options                                                        #
#############################################################################
case $ac_sys_system/$ac_sys_release in
Darwin/*)
    # Use -undefined dynamic_lookup whenever possible (10.3 and later).
    # This allows an extension to be used in any Python
    if test ${MACOSX_DEPLOYMENT_TARGET} '>' 10.2; then
        if test "${enable_universalsdk}"; then
            LDFLAGS="${UNIVERSAL_ARCH_FLAGS} -isysroot ${UNIVERSALSDK} ${LDFLAGS}"
        fi
    fi
    ;;
esac


#############################################################################
# Checks for header files                                                   #
#############################################################################
AC_HEADER_STDC


#############################################################################
# Checks for installer                                                      #
#############################################################################
AC_PROG_INSTALL


#############################################################################
# Checks for libtool.  We need to do this after the compiler has been       #
# determined.                                                               #
#############################################################################
AC_PROG_LIBTOOL


#############################################################################
# Checks for typedefs, structures, and compiler characteristics             #
#############################################################################
AC_TYPE_SIZE_T
# AC_TYPE_LONG_LONG_INT
AC_C_CONST
AC_C_INLINE


#############################################################################
# Checks for library functions                                              #
#############################################################################
# AC_CHECK_LIB(m, sincos, [AC_DEFINE([HAVE_SINCOS],[1],[Define to 1 if your system has `sincos'.])])


#############################################################################
# Checks for NaN/Inf checking option (G_NAN_CHECK)                          #
#############################################################################
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


#############################################################################
# Checks for Range checking option (G_RANGE_CHECK)                          #
#############################################################################
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


#############################################################################
# Checks for Small memory option (G_SMALL_MEMORY)                           #
#############################################################################
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


#############################################################################
# Checks if Python bindings are to be build                                 #
# ------------------------------------------------------------------------- #
# Python binding requires Python and either the files gammalib_wrap.cpp and #
# gammalib.py or swig to create these files.  The following conditionals    #
# are set:                                                                  #
# WITH_PYTHON - Python bindings will be built                               #
# HAVE_SWIG   - swig is available                                           #
#############################################################################
# Initialise flags
enable_python="no"
has_python="no"
has_python_h="no"
has_swig="no"

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

# Signal if Python is enabled
if [test "x$ac_enable_python_binding" = "xyes" -a "x$has_python" = "xyes" -a "x$has_python_h" = "xyes"]; then
  if [test "x$has_swig" = "xyes"]; then
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


#############################################################################
# Checks for readline and ncurses libraries                                 #
#############################################################################
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


#############################################################################
# Checks for cfitsio library                                                #
#############################################################################
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


#############################################################################
# Checks for optional instrument specific interfaces                        #
#############################################################################
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


#############################################################################
# Checks for Doxygen                                                        #
#############################################################################
DX_HTML_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_INIT_DOXYGEN([GammaLib], [doc/Doxyfile], [doc/doxygen])


#############################################################################
# Set list of output files to be created                                    #
#############################################################################
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


#############################################################################
# Print configuration summary to console                                    #
#############################################################################
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

# Dump Python bindings information
if test "x$ac_enable_python_binding" = "xyes"; then
  if test "x$has_python" = "xyes"; then
    if test "x$has_python_h" = "xyes"; then
      echo "  * Make Python bindings         (yes)"
    else
      echo "  - Make Python bindings         (no)    Python.h not installed"
    fi
  else
    echo "  - Make Python bindings         (no)    Python not installed"
  fi
else
  echo "  - Make Python bindings         (no)"
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
