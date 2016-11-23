.. _sec_xspec:

Xspec interface
---------------

Overview
~~~~~~~~

The following figure presents an overview over the C++ classes of the Xspec
module and their relations.

.. _fig_uml_xspec:

.. figure:: uml_xspec.png
   :width: 40%
   :align: center

   *Xspec module*

The Xspec module provides an interface to data in the XSPEC format. Event
data storage is implemented by the :doxy:`GPha` class that provides
the number of events measured per channel. The effective area taking
into account any event selection cuts is implemented by the :doxy:`GArf`
class. And energy redistribution information is implemented by the
:doxy:`GRmf` class.
