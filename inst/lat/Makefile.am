## Process this file with automake to produce Makefile.in

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
caldb_indx  = caldb/data/glast/lat/caldb.indx
caldb_add   = caldb/data/glast/lat/bcf/irf_index.fits
caldb_ea    = caldb/data/glast/lat/bcf/ea/aeff_P6_v3_diff_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P6_v3_diff_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7CLEAN_V6_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7CLEAN_V6_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7SOURCE_V6_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7SOURCE_V6_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7TRANSIENT_V6_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7TRANSIENT_V6_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7ULTRACLEAN_V6_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7ULTRACLEAN_V6_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_CLEAN_V15_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_CLEAN_V15_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_SOURCE_V15_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_SOURCE_V15_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_TRANSIENT_V15_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_TRANSIENT_V15_front.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_ULTRACLEAN_V15_back.fits \
              caldb/data/glast/lat/bcf/ea/aeff_P7REP_ULTRACLEAN_V15_front.fits
caldb_edisp = caldb/data/glast/lat/bcf/edisp/edisp_P6_v3_diff_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P6_v3_diff_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7CLEAN_V6_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7CLEAN_V6_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7SOURCE_V6_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7SOURCE_V6_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7TRANSIENT_V6_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7TRANSIENT_V6_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7ULTRACLEAN_V6_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7ULTRACLEAN_V6_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_CLEAN_V15_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_CLEAN_V15_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_SOURCE_V15_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_SOURCE_V15_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_TRANSIENT_V15_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_TRANSIENT_V15_front.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_ULTRACLEAN_V15_back.fits \
              caldb/data/glast/lat/bcf/edisp/edisp_P7REP_ULTRACLEAN_V15_front.fits
caldb_psf   = caldb/data/glast/lat/bcf/psf/psf_P6_v3_diff_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P6_v3_diff_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7CLEAN_V6_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7CLEAN_V6_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7SOURCE_V6_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7SOURCE_V6_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7TRANSIENT_V6_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7TRANSIENT_V6_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7ULTRACLEAN_V6_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7ULTRACLEAN_V6_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_CLEAN_V15_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_CLEAN_V15_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_SOURCE_V15_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_SOURCE_V15_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_TRANSIENT_V15_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_TRANSIENT_V15_front.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_ULTRACLEAN_V15_back.fits \
              caldb/data/glast/lat/bcf/psf/psf_P7REP_ULTRACLEAN_V15_front.fits

# Build libtool library
noinst_LTLIBRARIES = liblat.la
liblat_la_SOURCES = $(sources)
liblat_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST      = $(caldb_indx) $(caldb_add) $(caldb_ea) $(caldb_edisp) $(caldb_psf)
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
