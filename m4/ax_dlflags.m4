# ===========================================================================
#
# SYNOPSIS
#
#   AX_DLFLAGS
#
# DESCRIPTION
#
#   This macro determines the value of a DL flags such as RTLD_GLOBAL or
#   RTLD_NOW from the dlfcn.h file.
#
# USAGE
#
#   AX_DLFLAGS(RTLD_GLOBAL)
#   AC_MSG_RESULT($DLFLAGS_RTLD_GLOBAL)
#
# LICENSE
#
#   Copyright (c) 2021 Juergen Knoedlseder
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

#serial 8

AC_DEFUN([AX_DLFLAGS], [
  DLFLAGS_$1=""
  AC_CHECK_HEADERS([dlfcn.h])
  AS_IF([test "x$ac_cv_header_dlfcn_h" = "xyes"],[
    AC_MSG_CHECKING([for value of $1])
    for flags in 1 2 4 8 16 32 64 128 256 512 ; do
      AC_EGREP_CPP([test_true],[
        #include <dlfcn.h>
        #if defined($1) && $1 == $flags
        test_true
        #endif
        ],[DLFLAGS_$1=$flags],[])
    done
    ])
  AC_MSG_RESULT($DLFLAGS_$1)
  AC_SUBST([DLFLAGS_$1])
])
