## Process this file with automake to produce Makefile.in

# Include M4 macros
ACLOCAL_AMFLAGS = -I m4

# Detect configuration
if WITH_INST
INST_DIR = inst
endif

# Subdirectories
SUBDIRS = include $(INST_DIR) src test pyext

# Pkgconfig directory
pkgconfigdir = $(libdir)/pkgconfig

# Files to install in Pkgconfig directory
pkgconfig_DATA = gammalib.pc

# Manual pages to install
man_MANS = 

# Architecture independent data files installed in the package's data directory
doc_DATA = \
	   AUTHORS \
	   ChangeLog \
	   COPYING \
	   INSTALL \
	   LICENSE \
	   NEWS \
	   README

# Extra files to be distributed in the tar file
#EXTRA_DIST = autogen.sh \
#             gammalib.pc.in



# Automatically reconfigure libtool
LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck
