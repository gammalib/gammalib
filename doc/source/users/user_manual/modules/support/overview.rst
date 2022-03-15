Overview
~~~~~~~~

The following figure presents an overview over the classes of the support module
and their relations.

.. _fig_uml_support:

.. figure:: uml_support.png
   :width: 75%
   :align: center

   *Support module*


The support module contains a number of functions and classes that are 
needed to support core functionnalities of GammaLib. Most of the classes are
unrelated.

The :doxy:`GBilinear` class provides support for bilinear interpolation
in a table of values.

The :doxy:`GDaemon` class implements the GammaLib daemon that is started
in the background and that collects GammaLib usage statistics and information
about the carbon footprint of the software. Usage statistics is written by the
:doxy:`GApplication` class in a low level file situated at
``$HOME/.gamma/statistics.csv``. To avoid growing that this files becomes too
large the daemon creates high-level information from this file and writes them
in a more compact format to the file ``$HOME/.gamma/statistics.xml``.

The :doxy:`GCsv` class supports handling of column separated value tables.

The :doxy:`GExceptionHandler` class implements an exception handler that is
used through GammaLib for exception handling. The :doxy:`GException` class
derives from this handler, and implements a number of sub-classes that are
actually thrown in exceptions.

The :doxy:`GFilename` class supports handling of file names and FITS file
extensions. All file names in GammaLib are passed through instances of this
class.

The :doxy:`GNodeArray` implements methods for linear interpolation between node
values. This is the central interpolation class that is used in GammaLib.

The :doxy:`GRan` class implements a random number generator that is widely used
for Monte Carlo simulations.

The :doxy:`GTools.hpp` module is not a class, but a collection of constants and
functions that is widely used in GammaLib.

The abstract :doxy:`GUrl` base class represents a unified location for some
information that is independent of the media. The classes :doxy:`GUrlFile`
and :doxy:`GUrlString` are derived from :doxy:`GUrl` and implement a file storage 
and a string storage, respectively.
