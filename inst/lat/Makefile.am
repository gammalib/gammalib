#############################################################################
# Makefile.am for lat module                                                #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2018 Juergen Knoedlseder                               #
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

# Set sub directories
SUBDIRS = test

# Darwin support
MACOSX_DEPLOYMENT_TARGET=@CONFIGURE_MACOSX_DEPLOYMENT_TARGET@
@EXPORT_MACOSX_DEPLOYMENT_TARGET@export MACOSX_DEPLOYMENT_TARGET

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support \
           -I$(top_srcdir)/inst/lat/include

# Define sources for this directory
sources = src/GLATObservation.cpp \
          src/GLATEventList.cpp \
          src/GLATEventCube.cpp \
          src/GLATEventAtom.cpp \
          src/GLATEventBin.cpp \
          src/GLATResponse.cpp \
          src/GLATResponseTable.cpp \
          src/GLATAeff.cpp \
          src/GLATPsf.cpp \
          src/GLATPsfBase.cpp \
          src/GLATPsfV1.cpp \
          src/GLATPsfV3.cpp \
          src/GLATEdisp.cpp \
          src/GLATEfficiency.cpp \
          src/GLATMeanPsf.cpp \
          src/GLATInstDir.cpp \
          src/GLATRoi.cpp \
          src/GLATLtCube.cpp \
          src/GLATLtCubeMap.cpp \
          src/GLATException.cpp

# Define headers to be installed
pkginclude_HEADERS = include/GLATLib.hpp \
                     include/GLATObservation.hpp \
                     include/GLATEventList.hpp \
                     include/GLATEventCube.hpp \
                     include/GLATEventAtom.hpp \
                     include/GLATEventBin.hpp \
                     include/GLATResponse.hpp \
                     include/GLATResponseTable.hpp \
                     include/GLATAeff.hpp \
                     include/GLATPsf.hpp \
                     include/GLATPsfBase.hpp \
                     include/GLATPsfV1.hpp \
                     include/GLATPsfV3.hpp \
                     include/GLATEdisp.hpp \
                     include/GLATEfficiency.hpp \
                     include/GLATMeanPsf.hpp \
                     include/GLATInstDir.hpp \
                     include/GLATRoi.hpp \
                     include/GLATLtCube.hpp \
                     include/GLATLtCubeMap.hpp \
                     include/GLATException.hpp

# CALDB to distribute
caldb_indx  = caldb/data/glast/lat/caldb.indx.gz
caldb_add   = caldb/data/glast/lat/bcf/irf_index.fits.gz
caldb_ea    = caldb/data/glast/lat/bcf/ea/aeff_P6_v3_diff_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P6_v3_diff_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7CLEAN_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7CLEAN_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7SOURCE_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7SOURCE_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7TRANSIENT_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7TRANSIENT_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7ULTRACLEAN_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7ULTRACLEAN_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_CLEAN_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_CLEAN_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_SOURCE_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_SOURCE_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_TRANSIENT_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_TRANSIENT_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_ULTRACLEAN_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_ULTRACLEAN_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_CLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_CLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_CLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT010E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT010E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT010E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT010_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT010_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT010_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT015S_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT015S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT015S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT020E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT020E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT020E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT020_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT020_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT020_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100S_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_TRANSIENT100_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEANVETO_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEANVETO_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEANVETO_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_CLEAN_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_CLEAN_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_CLEAN_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCEVETO_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCEVETO_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCEVETO_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCE_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCE_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCE_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT010E_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT010E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT010E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT010_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT010_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT010_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT015S_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT015S_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT015S_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT020E_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT020E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT020E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT020_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT020_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT020_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100E_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100S_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100S_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100S_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_TRANSIENT100_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_ULTRACLEANVETO_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_ULTRACLEANVETO_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_ULTRACLEANVETO_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_ULTRACLEAN_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_ULTRACLEAN_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/ea/aeff_P8R3_ULTRACLEAN_V2_PSF.fits.gz
caldb_edisp = caldb/data/glast/lat/bcf/edisp/edisp_P6_v3_diff_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P6_v3_diff_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7CLEAN_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7CLEAN_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7SOURCE_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7SOURCE_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7TRANSIENT_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7TRANSIENT_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7ULTRACLEAN_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7ULTRACLEAN_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_CLEAN_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_CLEAN_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_SOURCE_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_SOURCE_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_TRANSIENT_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_TRANSIENT_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_ULTRACLEAN_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_ULTRACLEAN_V15_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_CLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_CLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_CLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010E_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT015S_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT015S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT015S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020E_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100E_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100S_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEANVETO_V6_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEANVETO_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEANVETO_V6_PSF.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_CLEAN_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_CLEAN_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_CLEAN_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCEVETO_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCEVETO_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCEVETO_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCE_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCE_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCE_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT010E_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT010E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT010E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT010_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT010_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT010_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT015S_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT015S_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT015S_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT020E_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT020E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT020E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT020_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT020_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT020_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100E_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100S_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100S_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100S_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_TRANSIENT100_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_ULTRACLEANVETO_V2_EDISP.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_ULTRACLEANVETO_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_ULTRACLEANVETO_V2_PSF.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_ULTRACLEAN_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_ULTRACLEAN_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R3_ULTRACLEAN_V2_PSF.fits.gz
caldb_psf   = caldb/data/glast/lat/bcf/psf/psf_P6_v3_diff_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P6_v3_diff_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7CLEAN_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7CLEAN_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7SOURCE_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7SOURCE_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7TRANSIENT_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7TRANSIENT_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7ULTRACLEAN_V6_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7ULTRACLEAN_V6_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_CLEAN_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_CLEAN_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_SOURCE_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_SOURCE_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_TRANSIENT_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_TRANSIENT_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_ULTRACLEAN_V15_back.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_ULTRACLEAN_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_CLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_CLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_CLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT010E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT010E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT010E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT010_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT010_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT010_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT015S_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT015S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT015S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT020E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT020E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT020E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT020_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT020_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT020_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100S_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_TRANSIENT100_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEANVETO_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEANVETO_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEANVETO_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_CLEAN_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_CLEAN_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_CLEAN_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCEVETO_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCEVETO_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCEVETO_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCE_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCE_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCE_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT010E_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT010E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT010E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT010_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT010_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT010_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT015S_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT015S_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT015S_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT020E_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT020E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT020E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT020_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT020_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT020_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100E_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100E_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100E_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100S_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100S_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100S_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_TRANSIENT100_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_ULTRACLEANVETO_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_ULTRACLEANVETO_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_ULTRACLEANVETO_V2_PSF.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_ULTRACLEAN_V2_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_ULTRACLEAN_V2_FB.fits.gz \
              caldb/data/glast/lat/bcf/psf/psf_P8R3_ULTRACLEAN_V2_PSF.fits.gz

# Distribute calibration database files
caldbindxdir    = $(datadir)/caldb/data/glast/lat
caldbindx_DATA  = $(caldb_indx)
caldbadddir     = $(datadir)/caldb/data/glast/lat/bcf
caldbadd_DATA   = $(caldb_add)
caldbeadir      = $(datadir)/caldb/data/glast/lat/bcf/ea
caldbea_DATA    = $(caldb_ea)
caldbedispdir   = $(datadir)/caldb/data/glast/lat/bcf/edisp
caldbedisp_DATA = $(caldb_edisp)
caldbpsfdir     = $(datadir)/caldb/data/glast/lat/bcf/psf
caldbpsf_DATA   = $(caldb_psf)

# Test data to distribute
testdatalat       = test/data/p8v2/binned_expmap.fits.gz \
                    test/data/p8v2/cntmap.fits.gz \
                    test/data/p8v2/ft1.fits.gz \
                    test/data/p8v2/ft2.fits.gz \
                    test/data/p8v2/iso_P8R2_SOURCE_V6_v06.txt \
                    test/data/p8v2/ltcube.fits.gz \
                    test/data/p8v2/modmap.fits.gz \
                    test/data/p8v2/obs_binned.xml \
                    test/data/p8v2/obs_unbinned.xml \
                    test/data/p8v2/source_model.xml \
                    test/data/p8v2/srcmap.fits.gz
testdatalatdir    = $(pyexecdir)/gammalib/tests/lat/data/p8v2
testdatalat_DATA  = $(testdatalat)
#
testcaldblat      = test/caldb/data/glast/lat/caldb.indx.gz
testcaldblatdir   = $(pyexecdir)/gammalib/tests/lat/caldb/data/glast/lat
testcaldblat_DATA = $(testcaldblat)
#
testbcflat        = test/caldb/data/glast/lat/bcf/irf_index.fits.gz \
                    test/caldb/data/glast/lat/bcf/valid_evclass_selections.txt	\
                    test/caldb/data/glast/lat/bcf/valid_evtype_selections.txt
testbcflatdir     = $(pyexecdir)/gammalib/tests/lat/caldb/data/glast/lat/bcf
testbcflat_DATA   = $(testbcflat)
#
testealat         = test/caldb/data/glast/lat/bcf/ea/aeff_P6_v3_diff_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P6_v3_diff_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P7SOURCE_V6_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P7SOURCE_V6_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P7REP_SOURCE_V15_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P7REP_SOURCE_V15_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_EDISP.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_FB.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P8R2_SOURCE_V6_PSF.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCE_V2_EDISP.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCE_V2_FB.fits.gz \
                    test/caldb/data/glast/lat/bcf/ea/aeff_P8R3_SOURCE_V2_PSF.fits.gz
testealatdir      = $(pyexecdir)/gammalib/tests/lat/caldb/data/glast/lat/bcf/ea
testealat_DATA    = $(testealat)
#
testedisplat      = test/caldb/data/glast/lat/bcf/edisp/edisp_P6_v3_diff_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P6_v3_diff_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P7SOURCE_V6_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P7SOURCE_V6_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P7REP_SOURCE_V15_back.fits \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P7REP_SOURCE_V15_front.fits \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_EDISP.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_FB.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_PSF.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCE_V2_EDISP.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCE_V2_FB.fits.gz \
                    test/caldb/data/glast/lat/bcf/edisp/edisp_P8R3_SOURCE_V2_PSF.fits.gz
testedisplatdir   = $(pyexecdir)/gammalib/tests/lat/caldb/data/glast/lat/bcf/edisp
testedisplat_DATA = $(testedisplat)
#
testpsflat        = test/caldb/data/glast/lat/bcf/psf/psf_P6_v3_diff_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P6_v3_diff_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P7SOURCE_V6_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P7SOURCE_V6_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P7REP_SOURCE_V15_back.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P7REP_SOURCE_V15_front.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_EDISP.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_FB.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P8R2_SOURCE_V6_PSF.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCE_V2_EDISP.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCE_V2_FB.fits.gz \
                    test/caldb/data/glast/lat/bcf/psf/psf_P8R3_SOURCE_V2_PSF.fits.gz
testpsflatdir     = $(pyexecdir)/gammalib/tests/lat/caldb/data/glast/lat/bcf/psf
testpsflat_DATA   = $(testpsflat)

# Build libtool library
noinst_LTLIBRARIES = liblat.la
liblat_la_SOURCES = $(sources)
liblat_la_LDFLAGS =

# Extra files to be distributed in the tar file
EXTRA_DIST = README.md caldb pyext \
             include/GLATLib.hpp \
             test/test_LAT.hpp test/test_LAT.py test/caldb test/data
