##### 
#
# SYNOPSIS
#
#   AC_FIND_LIB_HEADERS(LIBRARY, FUNCTION, HEADERS
#                       [ACTION-IF-FOUND],
#                       [ACTION-IF-NOT-FOUND],
#                       [OTHER-LIBRARIES])
#
# DESCRIPTION
#
#   This macro searches for a particular library on your system by
#   scanning a number of default directories. You can use the optional
#   third argument to search also for headers using AC_CHECK_HEADERS.
#
#   The directories that are searched are (in the given order):
#   $gammalib_prefix/lib64
#   $gammalib_prefix/lib
#   /opt/local/lib64       (for Mac OS X)
#   /opt/local/lib         (for Mac OS X)
#   /usr/local/lib64
#   /usr/local/lib
#   /usr/lib64
#   /usr/lib
#   /lib64                 (for readline/ncurses is some distros)
#   /lib                   (for readline/ncurses is some distros)
#
#   If the library is found, the function adds it's name to LIBS.
#   Furthermore, if the library is found in a specific directory, the
#   directory name is added to LDFLAGS. If the headers are found in a
#   specific directory, the directory name is added to CPPFLAGS.
#
#   In configure.in, use as:
#
#      AC_FIND_LIB_HEADERS([readline], [readline], [], [], [], [])
#
# LAST MODIFICATION
#
#   2011-07-13
#
# COPYLEFT
#
#   Copyright (c) 2011 Juergen Knodlseder <knodlseder@cesr.fr>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

AC_DEFUN([AC_FIND_LIB_HEADER],[
  #
  # Initialise status
  #
  AS_VAR_SET([ac_find_lib_status], [no])
  AS_VAR_SET([ac_find_header_status], [no])
  #
  # Signal checking
  #
  AC_MSG_CHECKING([for $2 in -l$1])
  #
  # Save original values
  #
  # echo "in: "$LIBS $LDFLAGS $CPPFLAGS
  ac_find_lib_save_LIBS=$LIBS
  ac_find_lib_save_LDFLAGS=$LDFLAGS
  ac_find_lib_add_LDFLAGS=
  ac_find_lib_base=
  #
  # Set link flags for testing
  #
  LIBS="-l$1 $6 $LIBS"
  #
  # Loop over all directories
  #
  for i in '' "$gammalib_prefix" /opt/local /usr/local /usr /; do
    #
    # Search library in lib64 first
    #
    if test -z "$i"; then
      LDFLAGS="$ac_find_lib_save_LDFLAGS"
      ac_message="found"
    else
      LDFLAGS="$ac_find_lib_save_LDFLAGS -L$i/lib64"
      ac_message="found in $i/lib64"
    fi
    AC_LINK_IFELSE([AC_LANG_CALL([], [$2])],
                 [AS_VAR_SET([ac_find_lib_status], [yes])],
                 [AS_VAR_SET([ac_find_lib_status], [no])])
    LDFLAGS=$ac_find_lib_save_LDFLAGS
    if test "x$ac_find_lib_status" = "xyes"; then
      if test -z "$i"; then
        ac_find_lib_add_LDFLAGS=
      else
        ac_find_lib_add_LDFLAGS=" -L$i/lib64"
        ac_find_lib_base="$i"
      fi
      break
    fi
    #
    # If still alive, check library in lib next
    #
    if test -z "$i"; then
      LDFLAGS="$ac_find_lib_save_LDFLAGS"
      ac_message="found"
    else
      LDFLAGS="$ac_find_lib_save_LDFLAGS -L$i/lib"
      ac_message="found in $i/lib"
    fi
    AC_LINK_IFELSE([AC_LANG_CALL([], [$2])],
                   [AS_VAR_SET([ac_find_lib_status], [yes])],
                   [AS_VAR_SET([ac_find_lib_status], [no])])
    LDFLAGS=$ac_find_lib_save_LDFLAGS
    if test "x$ac_find_lib_status" = "xyes"; then
      if test -z "$i"; then
        ac_find_lib_add_LDFLAGS=
      else
        ac_find_lib_add_LDFLAGS=" -L$i/lib"
        ac_find_lib_base="$i"
      fi
      break
    fi
  done
  #
  # Recover original values
  #
  LIBS=$ac_find_lib_save_LIBS
  LDFLAGS=$ac_find_lib_save_LDFLAGS
  #
  # Inform about test result
  #
  if test "x$ac_find_lib_status" = "xyes"; then
    AC_MSG_RESULT($ac_message)
  else
    AC_MSG_RESULT([not found])
  fi
  #
  # Optionally check for presence of headers (see AC_CHECK_HEADERS)
  #
  if test "x$3" != "x"; then
    ac_find_lib_save_CPPFLAGS=$CPPFLAGS
    ac_find_lib_add_CPPFLAGS=
    if test -z "$ac_find_lib_base"; then
      CPPFLAGS="$ac_find_lib_save_CPPFLAGS"
    else
      CPPFLAGS="$ac_find_lib_save_CPPFLAGS -I$ac_find_lib_base/include"
    fi
    AC_CHECK_HEADERS([$3],
                     [AS_VAR_SET([ac_find_header_status], [yes])
                      break])
    CPPFLAGS=$ac_find_lib_save_CPPFLAGS
    if test "x$ac_find_header_status" = "xyes"; then
      if test -z "$ac_find_lib_base"; then
        CPPFLAGS="$CPPFLAGS"
        AC_MSG_NOTICE([$1 headers found])
      else
        CPPFLAGS="$CPPFLAGS -I$ac_find_lib_base/include"
        AC_MSG_NOTICE([$1 headers found in $ac_find_lib_base/include])
      fi
    else
      AC_MSG_WARN([No $1 headers found])
    fi
    AS_IF([test "$ac_find_lib_status" != no -a "ac_find_header_status" != no ],
          [LIBS="-l$1 $LIBS"
           LDFLAGS="$LDFLAGS $ac_find_lib_add_LDFLAGS"
           $4], [$5])
  else
    AS_IF([test "$ac_find_lib_status" != no],
          [LIBS="-l$1 $LIBS"
           LDFLAGS="$LDFLAGS $ac_find_lib_add_LDFLAGS"
           $4], [$5])
  fi
  # echo "out: "$LIBS $LDFLAGS $CPPFLAGS
])
