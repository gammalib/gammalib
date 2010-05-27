## Process this file with automake to produce Makefile.in

# Include M4 macros
ACLOCAL_AMFLAGS = -I m4

# Subdirectories
SUBDIRS = include inst src test pyext

# Extra files to be distributed in the tar file
EXTRA_DIST = autogen.sh \
             gammalib.pc.in

# Install pkg-config information
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gammalib.pc

