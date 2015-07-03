
Overview
~~~~~~~~

The following figure presents an overview over the C++ classes of the CTA
module and their relations.

.. _fig_uml_cta:

.. figure:: uml_cta.png
   :width: 100%
   :align: center

   *CTA module*

The CTA module provides an instrument interface for Imaging Air Cherenkov 
Telescopes (IACT). The interface currently supports the event list data
format developed in the context of the Cherenkov Telescope Array (CTA)
project and provides support for various instrument response formats that
are currently used in the CTA consortium or that are under study. The interface
provides support for unbinned and binned maximum likelihood analysis, and
implements Monte Carlo simulations of high-level data based on the 
instrumental response functions.

Two types of observations are implemented so far:
:doxy:`GCTAObservation` that derives from :doxy:`GObservation` and that
describes either a binned or unbinned data set, and
:doxy:`GCTAOnOffObservation` that describes an On-Off observation
used in classical Cherenkov Telescope data analysis.
Binned or unbinned observations are collected in the
:doxy:`GObservations` container class, while the On-Off observations
are collected in a specific :doxy:`GCTAOnOffObservations` container
class (note that the :doxy:`GCTAOnOffObservations` container may vanish
in the future when :doxy:`GCTAOnOffObservation` objects will also be
collected by the :doxy:`GObservations` container class).

The data of binned observations are stored by the
:doxy:`GCTAEventCube` class. The :doxy:`GCTAEventCube` class holds the
binned event data in form of a sky map, implemented by the
:doxy:`GSkyMap` class. The sky coordinates of all sky map pixels
are stored in an array of CTA instrument directions, implemented
by the :doxy:`GCTAInstDir` class which holds a single :doxy:`GSkyDir`
object. The mean energies of the event cube
are stored in an array of :doxy:`GEnergy` objects, and the mean
time is stored by a :doxy:`GTime` object. The :doxy:`GCTAEventCube` class
holds in fact only a single event bin, implemented by the
:doxy:`GCTAEventBin` class. When the event bin is accessed using
the :doxy:`GCTAEventCube::operator[]` operator, the operator updates
references to the event cube data so that the :doxy:`GCTAEventBin`
object represents the selected event bin. This allows a memory-efficient
storage of event bin information (without storing for example the
instrument direction or the energy for each bin), while preserving
the abstract data model where an event cube of abstract type
:doxy:`GEventCube` is composed of event bins of abstract type 
:doxy:`GEventBin` (see :ref:`sec_obs`).

The data of unbinned observations are stored by the :doxy:`GCTAEventList`
class. The :doxy:`GCTAEventList` class is a container of :doxy:`GCTAEventAtom`
objects that represent individual events. Each event is composed of
a :doxy:`GCTAInstDir` object, a :doxy:`GEnergy` object and a :doxy:`GTime` object.
The region of interest covered by an event list is described by the
:doxy:`GCTARoi` class that derives from the abstract :doxy:`GRoi` class.

In addition to the event, :doxy:`GCTAObservation` holds pointing 
information implemented by the :doxy:`GCTAPointing` class.
:doxy:`GCTAObservation` holds furthermore response information.
Response information derives from the abstract :doxy:`GCTAResponse` class
that is impleted either as a factorized Instrument Response
Function by :doxy:`GCTAResponseIrf` or as a precomputed response for
stacked analysis by :doxy:`GCTAResponseCube`.
For stacked analysis the response is composed of the exposure
(implemented by :doxy:`GCTACubeExposure`), the livetime averaged point
spread function (implemented by :doxy:`GCTACubePsf`) and the livetime
averaged background rate (implemented by :doxy:`GCTACubeBackground`).
For all other analysis, the factorized IRFs are used directly.
Dependent on the response format, the components of the factorization
are stored in classes that derive from the abstract :doxy:`GCTAAeff`, 
:doxy:`GCTAPsf`, :doxy:`GCTAEdisp`, and :doxy:`GCTABackground` classes
(see :ref:`sec_cta_response` for more details on the response
implementation).

Different models of the instrumental background are provided by the
:doxy:`GCTAModelIrfBackground`, :doxy:`GCTAModelCubeBackground`, and
:doxy:`GCTAModelRadialAcceptance` classes.
All these classes implement the abstract :doxy:`GModelData` base
class.
:doxy:`GCTAModelIrfBackground` and :doxy:`GCTAModelCubeBackground`
will use the background rate information that is attached to the 
response function for background modelling.
Specifically, :doxy:`GCTAModelIrfBackground` will use information
provided by :doxy:`GCTABackground` and :doxy:`GCTAModelCubeBackground`
will use information provided by :doxy:`GCTACubeBackground`
(see :ref:`sec_cta_background` for details on the background
model implementation).
