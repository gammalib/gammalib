## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

# Define sources for this directory
sources = GVector.cpp \
	  GMatrixTools.cpp \
	  GMatrixBase.cpp \
	  GMatrix.cpp \
	  GSymMatrix.cpp \
	  GSparseMatrix.cpp \
	  GSparseSymbolic.cpp \
	  GSparseNumeric.cpp

# Build libtool library
noinst_LTLIBRARIES = liblinalg.la
liblinalg_la_SOURCES = $(sources)
liblinalg_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
