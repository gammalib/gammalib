## Process this file with automake to produce Makefile.in

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src/support

# Define sources for this directory
sources = GXml.cpp \
          GXmlNode.cpp \
          GXmlDocument.cpp \
          GXmlText.cpp \
          GXmlElement.cpp \
          GXmlComment.cpp \
          GXmlAttribute.cpp \
          GException_xml.cpp

# Build libtool library
noinst_LTLIBRARIES = libxml.la
libxml_la_SOURCES = $(sources)
libxml_la_LDFLAGS = -version-info $(GAMMALIB_LT_VERSION)
