## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_builddir)/src \
	   -I$(top_srcdir)/src

# Define sources for this directory
sources = GOptimizer.cpp \
	  GOptimizerLM.cpp \
	  GOptimizerPars.cpp \
	  GOptimizerFunction.cpp
	
# Build libtool library
noinst_LTLIBRARIES = libopt.la
libopt_la_SOURCES = $(sources)
libopt_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
