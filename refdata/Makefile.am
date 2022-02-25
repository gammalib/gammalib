#############################################################################
# Makefile for reference data                                               #
# ------------------------------------------------------------------------- #
# Copyright (C) 2022 Juergen Knoedlseder                                    #
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

# Models to install
refdata = comptel_db_tofIV.txt \
          emission-factors.fits \
          ephem_jpl_de200.fits

# Install models
refdatadir   = $(datadir)/refdata
refdata_DATA = $(refdata)

# Extra files to be distributed in the tar file
EXTRA_DIST = $(refdata)

