## Process this file with automake to produce Makefile.in

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
              caldb/data/glast/lat/bcf/ea/aeff_P8R2_ULTRACLEAN_V6_PSF.fits.gz
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
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_ULTRACLEAN_V15_front.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_CLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_CLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_CLEAN_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_SOURCE_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT010_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT015S_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT015S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT015S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT020_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100E_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100E_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100E_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100S_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100S_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100S_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_TRANSIENT100_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEANVETO_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEANVETO_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEANVETO_V6_PSF.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEAN_V6_EDISP.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEAN_V6_FB.fits.gz \
              caldb/data/glast/lat/bcf/edisp/edisp_P8R2_ULTRACLEAN_V6_PSF.fits.gz
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
              caldb/data/glast/lat/bcf/psf/psf_P8R2_ULTRACLEAN_V6_PSF.fits.gz

# Build libtool library
noinst_LTLIBRARIES = liblat.la
liblat_la_SOURCES = $(sources)
liblat_la_LDFLAGS =

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

# Extra files to be distributed in the tar file
EXTRA_DIST = README.md caldb pyext \
             include/GLATLib.hpp \
             test/test_LAT.hpp test/test_LAT.py test/data
