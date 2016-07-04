#############################################################################
# Makefile.am for inst module                                               #
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

# Detect configuration
if WITH_INST_CTA
INST_CTA = cta
endif
if WITH_INST_LAT
INST_LAT = lat
endif
if WITH_INST_MWL
INST_MWL = mwl
endif
if WITH_INST_COM
INST_COM = com
endif

# Set sub directories
SUBDIRS = $(INST_MWL) $(INST_CTA) $(INST_LAT) $(INST_COM)

# Extra files to be distributed in the tar file
EXTRA_DIST = README.md
