.. _sec_vo:

Virtual Observatory interface
-----------------------------

Overview
~~~~~~~~

The VO module provides classes that interface GammaLib to the Virtual 
Observatory infrastructure. So far the module contains the :doxy:`GVOClient`
and the :doxy:`GVOHub` classes. :doxy:`GVOClient` implements a 
Virtual Observatory client while :doxy:`GVOHub` implement a hub for
inter-client communication.
If a client needs to connect to a Hub and no Hub is currently running
on the system, the client will launch it's own Hub.

:ref:`fig_uml_vo` presents an overview over the C++ classes of the VO
module and their relations.

.. _fig_uml_vo:

.. figure:: uml_vo.png
   :width: 30%

   VO module
