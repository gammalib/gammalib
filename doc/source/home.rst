Home
====

GammaLib is a self-contained, instrument independent, open source,
multi-platform C++ library that implements all code required for
high-level science analysis of astronomical gamma-ray data. GammaLib is
also wrapped into a Python module.

GammaLib does not rely on any third-party software, except of HEASARC's
cfitsio library that is used to implement the FITS interface. Large
parts of the code treat gamma-ray observations in an abstract
representation, and do neither depend on the characteristics of the
employed instrument, nor on the particular formats in which data and
instrument response functions are delivered. Instrument specific
aspects are implemented as isolated and well defined modules that
interact with the rest of the library through a common interface. This
philosophy also enables the joint analysis of data from different
instruments, providing a framework that allows for consistent
broad-band spectral fitting or imaging. So far, GammaLib supports
analysis of COMPTEL, Fermi-LAT, and Cherenkov telescope data (H.E.S.S., MAGIC, CTA).

GammaLib is free software distributed under the
`GNU GPL license version 3 <http://www.gnu.org/licenses/gpl.html>`_