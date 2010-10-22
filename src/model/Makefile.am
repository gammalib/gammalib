## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

# Define sources for this directory
sources = GModels.cpp \
	  GModel.cpp \
	  GModelPar.cpp \
	  GModelSpatial.cpp \
	  GModelSpatialPtsrc.cpp \
	  GModelSpectral.cpp \
	  GModelSpectralPlaw.cpp \
	  GModelTemporal.cpp \
	  GModelTemporalConst.cpp \
          GException_model.cpp

# Build libtool library
noinst_LTLIBRARIES = libmodel.la
libmodel_la_SOURCES = $(sources)
libmodel_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
