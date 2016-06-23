## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

# Program name
bin_PROGRAMS = createfits

# Program definition
createfits_SOURCES = createfits.cpp
createfits_LDFLAGS = @LDFLAGS@
createfits_CPPFLAGS = @CPPFLAGS@
createfits_LDADD = $(top_srcdir)/src/libgamma.la

