# GammaLib docs

This directory contains the gammalib documentation.

It consists of two parts:

1. Hand-written documentation in RST (restructured text) format
   that can be converted to HTML or PDF with a tool called
   Sphinx --- http://sphinx-doc.org
2. Automatic API (application programming interface) documentation
   that is contained in the C++ source code (.hpp and .cpp files)
   and can be extracted as HTML or PDF with a tool called 
   Doxygen --- http://doxygen.org

To generate the documentation in html format locally use these commands
in the top-level folder of the GammaLib repository:

	make doc

or if you only want the Sphinx or Doxygen part use:

	make sphinx
	make doxygen

Then point your web browser to

	doc/html/index.html
	doc/html/doxygen/index.html
