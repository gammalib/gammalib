#############################################################################
# Makefile.am for cta/test module                                           #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016 Juergen Knoedlseder                                    #
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
#############################################################################

# Test sources
test_CTA_SOURCES = test_CTA.cpp test_CTA.hpp
test_CTA_LDFLAGS = @LDFLAGS@
test_CTA_CPPFLAGS = @CPPFLAGS@ -I$(srcdir)/../include -I$(top_srcdir)/include
test_CTA_LDADD = $(top_builddir)/src/libgamma.la

# Test programs (list of executables that will be compiled)
check_PROGRAMS = test_CTA
