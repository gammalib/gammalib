.. _sec_app:

Creation of ftools applications
-------------------------------

Overview
~~~~~~~~

The following figure presents an overview over the C++ classes of the 
application module and their relations.

.. _fig_uml_app:

.. figure:: uml_app.png
   :width: 35%
   :align: center

   *Application module*

The application module provides classes to create ftools applications
using GammaLib. A ftool application provides a parameter interface
following the IRAF standard that is widely used in high-energy
astrophysics. The application is represented by the abstract
:doxy:`GApplication` base class. It contains a parameter container
implemented by the :doxy:`GApplicationPars` class that contains application
parameter implemented by the :doxy:`GApplicationPar` class. In addition, an
application contains a logger, implemented by the :doxy:`GLog` class, which 
allows message passing to the console and/or a dedicated log file.

