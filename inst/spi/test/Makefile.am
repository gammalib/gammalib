#############################################################################
# Makefile.am for spi/test module                                           #
# ------------------------------------------------------------------------- #
# Copyright (C) 2020 by Juergen Knoedlseder                                 #
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
test_SPI_SOURCES = test_SPI.cpp test_SPI.hpp
test_SPI_LDFLAGS = @LDFLAGS@
test_SPI_CPPFLAGS = @CPPFLAGS@ -I$(srcdir)/../include -I$(top_srcdir)/include
test_SPI_LDADD = $(top_builddir)/src/libgamma.la

# Test programs (list of executables that will be compiled)
check_PROGRAMS = test_SPI
