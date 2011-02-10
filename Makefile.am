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

#
# Add Doxygen rules
#
if DX_COND_doc

# HTML
if DX_COND_html
DX_CLEAN_HTML = @DX_DOCDIR@/html
endif DX_COND_html

# CHI & CHI
if DX_COND_chm
DX_CLEAN_CHM = @DX_DOCDIR@/chm
if DX_COND_chi
DX_CLEAN_CHI = @DX_DOCDIR@/@PACKAGE@.chi
endif DX_COND_chi
endif DX_COND_chm

# MAN
if DX_COND_man
DX_CLEAN_MAN = @DX_DOCDIR@/man
endif DX_COND_man

# RTF
if DX_COND_rtf
DX_CLEAN_RTF = @DX_DOCDIR@/rtf
endif DX_COND_rtf

# XML
if DX_COND_xml
DX_CLEAN_XML = @DX_DOCDIR@/xml
endif DX_COND_xml

# PS
if DX_COND_ps
DX_CLEAN_PS = @DX_DOCDIR@/@PACKAGE@.ps
DX_PS_GOAL = doxygen-ps
doxygen-ps: @DX_DOCDIR@/@PACKAGE@.ps
@DX_DOCDIR@/@PACKAGE@.ps: @DX_DOCDIR@/@PACKAGE@.tag
	   cd @DX_DOCDIR@/latex; \
	   rm -f *.aux *.toc *.idx *.ind *.ilg *.log *.out; \
	   $(DX_LATEX) refman.tex; \
	   $(MAKEINDEX_PATH) refman.idx; \
	   $(DX_LATEX) refman.tex; \
	   countdown=5; \
	   while $(DX_EGREP) 'Rerun (LaTeX|to get cross-references right)' \
	                     refman.log > /dev/null 2>&1 \
	      && test $$countdown -gt 0; do \
	       $(DX_LATEX) refman.tex; \
	       countdown=`expr $$countdown - 1`; \
	   done; \
	   $(DX_DVIPS) -o ../@PACKAGE@.ps refman.dvi
endif DX_COND_ps

# PDF
if DX_COND_pdf
DX_CLEAN_PDF = @DX_DOCDIR@/@PACKAGE@.pdf
DX_PDF_GOAL = doxygen-pdf
doxygen-pdf: @DX_DOCDIR@/@PACKAGE@.pdf
@DX_DOCDIR@/@PACKAGE@.pdf: @DX_DOCDIR@/@PACKAGE@.tag
	   cd @DX_DOCDIR@/latex; \
	   rm -f *.aux *.toc *.idx *.ind *.ilg *.log *.out; \
	   $(DX_PDFLATEX) refman.tex; \
	   $(DX_MAKEINDEX) refman.idx; \
	   $(DX_PDFLATEX) refman.tex; \
	   countdown=5; \
	   while $(DX_EGREP) 'Rerun (LaTeX|to get cross-references right)' \
	                     refman.log > /dev/null 2>&1 \
	      && test $$countdown -gt 0; do \
	       $(DX_PDFLATEX) refman.tex; \
	       countdown=`expr $$countdown - 1`; \
	   done; \
	   mv refman.pdf ../@PACKAGE@.pdf
endif DX_COND_pdf

# LaTeX
if DX_COND_latex
DX_CLEAN_LATEX = @DX_DOCDIR@/latex
endif DX_COND_latex

# Rules
.PHONY: doxygen-run doxygen-doc $(DX_PS_GOAL) $(DX_PDF_GOAL)

.INTERMEDIATE: doxygen-run $(DX_PS_GOAL) $(DX_PDF_GOAL)

doxygen-run: @DX_DOCDIR@/@PACKAGE@.tag

# Make Doxygen files
doxygen: doxygen-run $(DX_PS_GOAL) $(DX_PDF_GOAL)

# Install Doxygen files
DOCDIRS = doxygen/html doxygen/man/man3
doxygen-install:
	$(mkinstalldirs) $(DESTDIR)$(docdir)
	@for dir in $(DOCDIRS) ; do \
	  $(mkinstalldirs) $(DESTDIR)$(docdir)/$$dir && \
	  for file in $(srcdir)/doc/$$dir/*.* ; do \
	    echo "$(INSTALL_DATA) $$file $(DESTDIR)$(docdir)/$$dir" ; \
	    $(INSTALL_DATA) "$$file" $(DESTDIR)$(docdir)/$$dir ; \
	  done ; \
	done

# Uninstall Doxygen files
doxygen-uninstall:
	for dir in $(DOCDIRS) ; do \
	  rm -f -r $(DESTDIR)$(docdir)/$$dir ; \
	done

@DX_DOCDIR@/@PACKAGE@.tag: $(DX_CONFIG) $(pkginclude_HEADERS)
	    rm -rf @DX_DOCDIR@
	    $(DX_ENV) $(DX_DOXYGEN) $(srcdir)/$(DX_CONFIG)

DX_CLEANFILES = \
     @DX_DOCDIR@/@PACKAGE@.tag \
     -r \
     $(DX_CLEAN_HTML) \
     $(DX_CLEAN_CHM) \
     $(DX_CLEAN_CHI) \
     $(DX_CLEAN_MAN) \
     $(DX_CLEAN_RTF) \
     $(DX_CLEAN_XML) \
     $(DX_CLEAN_PS) \
     $(DX_CLEAN_PDF) \
     $(DX_CLEAN_LATEX)

endif DX_COND_doc
