#
# SYNOPSIS
#
#   AC_COUNTRY_CODE()
#
# DESCRIPTION
#
#   Note: Returns country code using curl command.
#
# LAST MODIFICATION
#
#   2022-02-21
#
# COPYLEFT
#
#   Copyright (c) 2022 Juergen Knodlseder <jknodlseder@irap.omp.eu>
#
#   This program is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see
#   <http://www.gnu.org/licenses/>.

AC_DEFUN([AC_COUNTRY_CODE],[
    COUNTRY_CODE=""
    AC_PATH_PROG([CURL],[curl])
    AC_PATH_PROG([TIMEOUT],[timeout])
    AC_PATH_PROG([PERL],[perl])
    if test -n "$CURL" ; then
        AC_MSG_CHECKING(country code)
        COUNTRY_CODE=`$CURL --silent http://ip-api.com/line/?fields=countryCode 2>/dev/null`
        AC_MSG_RESULT([$COUNTRY_CODE])
        AC_SUBST(COUNTRY_CODE)
        AC_DEFINE_UNQUOTED([G_COUNTRY_CODE], "${COUNTRY_CODE}", [Define country code])
        AC_DEFINE([G_HAS_CURL], [1], [curl tool available])
        AC_DEFINE_UNQUOTED([G_CURL], "${CURL}", [curl tool])
        
    fi
    if test -n "$TIMEOUT" ; then
        AC_DEFINE([G_HAS_TIMEOUT], [1], [timeout tool available])
        AC_DEFINE_UNQUOTED([G_TIMEOUT], "${TIMEOUT}", [timeout tool])
    fi
    if test -n "$PERL" ; then
        AC_DEFINE([G_HAS_PERL], [1], [perl available])
        AC_DEFINE_UNQUOTED([G_PERL], "${PERL}", [perl tool])
    fi
])
