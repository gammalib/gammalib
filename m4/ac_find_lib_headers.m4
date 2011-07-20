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

m4_define(AC_DEFUN([AC_FIND_LIB_HEADER],[
  # Initialise status
  AS_VAR_SET([ac_find_lib_status], [no])
  AS_VAR_SET([ac_find_header_status], [no])

  # Signal checking
  AC_MSG_CHECKING([for $2 in -l$1])

  # Save original values
  # echo "in: "$LIBS $LDFLAGS $CPPFLAGS
  ac_find_lib_save_LIBS=$LIBS
  ac_find_lib_save_LDFLAGS=$LDFLAGS
  ac_find_lib_save_CPPFLAGS=$CPPFLAGS
  ac_find_lib_add_LDFLAGS=
  ac_find_lib_add_CPPFLAGS=

  # Set link flags for testing
  LIBS="-l$1 $6 $LIBS"

  # Loop over all directories
  for i in '' "$gammalib_prefix" /opt/local /usr/local /usr; do
    for l in lib64 lib; do

      # Case A: headers are required
      if test "x$3" != "x"; then
        for header in $3; do
          if test "x$i" == "x"; then
            ac_find_lib_add_LDFLAGS=
            ac_find_lib_add_CPPFLAGS=
            LDFLAGS="$ac_find_lib_save_LDFLAGS"
            CPPFLAGS="$ac_find_lib_save_CPPFLAGS"
            ac_message="found"
          else
            ac_find_lib_add_LDFLAGS=" -L$i/$l"
            ac_find_lib_add_CPPFLAGS=" -I$i/include"
            LDFLAGS="$ac_find_lib_save_LDFLAGS$ac_find_lib_add_LDFLAGS"
            CPPFLAGS="$ac_find_lib_save_CPPFLAGS$ac_find_lib_add_CPPFLAGS"
            ac_message="found in $i/$lib and $i/$include"
          fi

          # Test if we can link the library
          AC_LINK_IFELSE([AC_LANG_CALL([], [$2])],
                         [AS_VAR_SET([ac_find_lib_status], [yes])],
                         [AS_VAR_SET([ac_find_lib_status], [no])])

          # Test if the preprocessor can handle the header
          AC_PREPROC_IFELSE([AC_LANG_SOURCE([@%:@include <$header>])],
                            [AS_VAR_SET([ac_find_header_status], [yes])],
                            [AS_VAR_SET([ac_find_header_status], [no])])

          # Test if we can compile the header (not used as problems with
          # readline
          #AC_COMPILE_IFELSE([AC_LANG_SOURCE([@%:@include <$header>])],
          #                  [AS_VAR_SET([ac_find_header_status], [yes])],
          #                  [AS_VAR_SET([ac_find_header_status], [no])])

          # Recover old flags
          LDFLAGS=$ac_find_lib_save_LDFLAGS
          CPPFLAGS=$ac_find_lib_save_CPPFLAGS
          if test "x$ac_find_lib_status" = "xyes" -a "x$ac_find_header_status" = "xyes"; then
            break
          else
            ac_find_lib_status="no"
          fi
        done
        if test "x$ac_find_lib_status" = "xyes"; then
          break
        fi
     
      # Case B: no headers are required
      else
        if test -z "$i"; then
          ac_find_lib_add_LDFLAGS=
          LDFLAGS="$ac_find_lib_save_LDFLAGS"
          ac_message="found"
        else
          ac_find_lib_add_LDFLAGS=" -L$i/$l"
            LDFLAGS="$ac_find_lib_save_LDFLAGS$ac_find_lib_add_LDFLAGS"
          ac_message="found in $ac_find_lib_add_LDFLAGS"
        fi
        AC_LINK_IFELSE([AC_LANG_CALL([], [$2])],
                       [AS_VAR_SET([ac_find_lib_status], [yes])],
                       [AS_VAR_SET([ac_find_lib_status], [no])])
        LDFLAGS=$ac_find_lib_save_LDFLAGS
        if test "x$ac_find_lib_status" = "xyes"; then
          break
        fi
      fi

    done
    if test "x$ac_find_lib_status" = "xyes"; then
      break
    fi
  done

  # Recover original LIBS value
  LIBS=$ac_find_lib_save_LIBS
  
  # Set flags
  if test "x$3" != "x"; then
    AS_IF([test "$ac_find_lib_status" == yes ],
          [LIBS="-l$1 $LIBS"
           LDFLAGS="$LDFLAGS$ac_find_lib_add_LDFLAGS"
           CPPFLAGS="$CPPFLAGS$ac_find_lib_add_CPPFLAGS"
           AC_MSG_RESULT($ac_message)
           $4],
          [AC_MSG_RESULT([not found])
           $5])
  else
    AS_IF([test "$ac_find_lib_status" == yes ],
          [LIBS="-l$1 $LIBS"
           LDFLAGS="$LDFLAGS$ac_find_lib_add_LDFLAGS"
           AC_MSG_RESULT($ac_message)
           $4],
          [AC_MSG_RESULT([not found])
           $5])
  fi
  # echo "out: "$LIBS $LDFLAGS $CPPFLAGS
]))
