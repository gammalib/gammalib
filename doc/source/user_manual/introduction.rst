.. _um_intro:

Introduction
============

Scope
-----

This User Manual provides a complete description of the GammaLib C++ library
and its Python module. It is equally well suited for beginners as for
experienced users.

Overview
--------

GammaLib is a self-contained, instrument independent, open source,
multi-platform C++ library and Python module that implements all code
required for the scientific analysis of astronomical gamma-ray data.
GammaLib works with high-level data, which are data that are already
calibrated and that generally come in form of event lists or histograms.

Except for HEASARC's cfitsio library that is used to implement the FITS
interface, GammaLib does not rely on any other third-party software.
This makes GammaLib basically independent of any other software
package, increasing the maintainability and enhancing the portability of
the library.

GammaLib potentially supports any gamma-ray astronomy instrument. Large
parts of the code treat gamma-ray observations in an abstract representation,
and do neither depend on the characteristics of the employed instrument,
nor on the particular formats in which data and instrument response
functions are delivered. Instrument specific aspects are implemented as
isolated and well defined modules that interact with the rest of the
library through a common interface. This philosophy also enables the
joint analysis of data from different instruments, providing a framework
that allows for consistent broad-band spectral fitting or imaging.

GammaLib source code is freely available under the GNU General Public
license version 3. Relevant links for download, issues and documentation
are:

* Download: http://cta.irap.omp.eu/ctools/download.html
* Issue tracker: https://cta-redmine.irap.omp.eu/projects/gammalib/issues
* Documentation: http://cta.irap.omp.eu/gammalib/

The present document applies to GammaLib version 0.11.

GammaLib is designed to compile on any POSIX compliant platform. So far, 
GammaLib has been successfully compiled and tested on Mac OS X, OpenBSD, OpenSolaris
(using the gcc compiler) and many Linux flavours. Pre-packed binary
versions of the code are also available for Mac OS X. For known problems
with specific platforms, please refer to the :ref:`issues`
section.

GammaLib makes heavily use of C++ classes. Instrument independency is achieved
by using abstract virtual base classes, which are implemented as derived
classes in the instrument specific modules.

GammaLib is organized into four software layers, each of which comprises a
number of modules (see :ref:`fig_structure`; the quoted module names
correspond to the folders in the source code distribution):

-  **High-level analysis support**

   GammaLib implements classes needed for the instrument independent high-level
   analysis of gamma-ray data, enabling the joint multi-instrument
   spectral and spatial fitting by forwards folding of parametric source
   and background models. This layer comprises modules for instrument
   independent handling of observations (obs), sky and background models
   (model), sky maps and sky coordinates (sky), and for the generation
   of analysis executables in form of ftools (app).

-  **Core services**

   GammaLib comprise modules for numerical computations (numerics), linear
   algebra (linalg), parameter optimization (opt), and support functions
   and classes (support).

-  **Interfaces**

   are implemented for reading and writing of FITS files (fits), XML
   files (xml), and Virtual Observatory interoperability (vo).

-  **Instrument specific modules**

   GammaLib support the analysis of data from the Cherenkov Telescope Array
   (cta), the Fermi/LAT telescope (lat), the COMPTEL telescope (com),
   and any multi-wavelength information in form of spectral data points (mwl).

.. _fig_structure:

.. figure:: structure.png
   :width: 75%
   :align: center

   *GammaLib structure*

GammaLib is developed by a team of enthousiastic gamma-ray 
astronomers with support from engineers.
We regularily organise
`coding sprints <https://cta-redmine.irap.omp.eu/projects/ctools/wiki/Coding_sprints>`_
where key developers but also newcomers meet to discuss the developments 
and next steps, and advance with the coding of the software.

The development of GammaLib has been initiated by scientists from `IRAP (Institut
de Recherche en Astrophysique et Planetologie) <http://www.irap.omp.eu/>`_, an
astrophysics laboratory of CNRS and of the `University Paul Sabatier 
<http://www.univ-tlse3.fr/>`_ situated in
Toulouse, France. GammaLib is based on past experience gained in developing
software for gamma-ray space missions, such as the COMPTEL telescope
aboard CGRO, the SPI telescope aboard INTEGRAL, and the LAT
telescope aboard Fermi. Today, the development of GammaLib is mainly driven
by the needs in ground-based gamma-ray astronomy, and in particular by the
development of the CTA observatory.

