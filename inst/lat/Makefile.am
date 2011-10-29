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
          src/GLATEdisp.cpp \
          src/GLATMeanPsf.cpp \
          src/GLATInstDir.cpp \
          src/GLATRoi.cpp \
          src/GLATPointing.cpp \
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
                     include/GLATEdisp.hpp \
                     include/GLATMeanPsf.hpp \
                     include/GLATInstDir.hpp \
                     include/GLATRoi.hpp \
                     include/GLATPointing.hpp \
                     include/GLATLtCube.hpp \
                     include/GLATLtCubeMap.hpp \
                     include/GLATException.hpp

# Calibration files to distribute
irfs = caldb/aeff_P6_v3_diff_front.fits \
       caldb/aeff_P6_v3_diff_back.fits \
       caldb/edisp_P6_v3_diff_front.fits \
       caldb/edisp_P6_v3_diff_back.fits \
       caldb/psf_P6_v3_diff_front.fits \
       caldb/psf_P6_v3_diff_back.fits

# CALDB to tdistribute
caldb_indx = caldb/data/fermi/lat/caldb.indx
caldb_ea = caldb/data/fermi/lat/bcf/ea/aeff_P6_v3_diff_front.fits \
           caldb/data/fermi/lat/bcf/ea/aeff_P6_v3_diff_back.fits
caldb_edisp = caldb/data/fermi/lat/bcf/edisp/edisp_P6_v3_diff_front.fits \
              caldb/data/fermi/lat/bcf/edisp/edisp_P6_v3_diff_back.fits
caldb_psf = caldb/data/fermi/lat/bcf/psf/psf_P6_v3_diff_front.fits \
            caldb/data/fermi/lat/bcf/psf/psf_P6_v3_diff_back.fits

# Build libtool library
noinst_LTLIBRARIES = liblat.la
liblat_la_SOURCES = $(sources)
liblat_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)

# Distribute calibration database files
EXTRA_DIST = $(irfs) $(caldb_indx) $(caldb_ea) $(caldb_edisp) $(caldb_psf)
caldbdir = $(datadir)/caldb/lat
caldb_DATA = $(irfs)
caldbindxdir = $(datadir)/caldb/data/fermi/lat
caldbindx_DATA = $(caldb_indx)
caldbeadir = $(datadir)/caldb/data/fermi/lat/bcf/ea
caldbea_DATA = $(caldb_ea)
caldbedispdir = $(datadir)/caldb/data/fermi/lat/bcf/edisp
caldbedisp_DATA = $(caldb_edisp)
caldbpsfdir = $(datadir)/caldb/data/fermi/lat/bcf/psf
caldbpsf_DATA = $(caldb_psf)
