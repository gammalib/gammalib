Overview
========

The following figure presents an overview over the classes of the obs
module and their relations.

.. _fig_uml_obs:

.. figure:: uml_obs.png
   :align: center
   :width: 100%

   *Overview over the obs module*

The central class of the obs module is the abstract base class
:doxy:`GObservation` which defines the instrument-independent interface for a
gamma-ray observation. A gamma-ray observation is defined for a single
specific instrument, and describes a time period during which the
instrument is in a given stable configuration that can be characterized
by a single specific response function. Each gamma-ray observation is
composed of events and a response function.

Observations are collected in the container class :doxy:`GObservations`
which is composed of a list of :doxy:`GObservation` elements (the list is of
arbitrary length; an empty list is a valid state of the :doxy:`GObservations`
class). The observation container is furthermore composed of a :doxy:`GModels`
model container class that holds a list of models used to describe the
event distributions of the observations (see :ref:`um_model`). The
:doxy:`GObservations` class presents the central element of all scientific data
analyses, as it combines all data and all models in a single entity.

Instrument specific implementations of :doxy:`GObservation` objects are
registered in the registry class :doxy:`GObservationRegistry` which
statically collects one instance of each instrument-specific observation
class that is available in GammaLib (see :ref:`um_registry` for a general
description of registry classes).

The instrument response for a given observation is defined by the
abstract base class :doxy:`GResponse`. This class is composed of the class
:doxy:`GCaldb` which implements the calibration data base that is required to
compute the response function for a given instrument and observation.
:doxy:`GCaldb` supports the HEASARC CALDB format
(http://heasarc.nasa.gov/docs/heasarc/caldb/), but is sufficiently
general to support also other formats (see :ref:`um_obs_caldb` to learn
how to setup and to use a calibration database).
In addition, the class contains members of :doxy:`GResponseCache` and
:doxy:`GResponseVectorCache` which implement value caching that avoids recomputation
of response information.

The events for a given observation are defined by the abstract base
class :doxy:`GEvents`. This class is composed of the classes :doxy:`GGti` and
:doxy:`GEbounds`. :doxy:`GGti` implements so called *Good Time Intervals*, which defines
the time period(s) during which the data were taken (see :ref:`um_obs_time`).
:doxy:`GEbounds` implements so called *Energy Boundaries*, which
define the energy intervals that are covered by the data (see 
:ref:`um_obs_energy`).

:doxy:`GEvents` is also a container for the individual events, implemented by the
abstract :doxy:`GEvent` base class. 
GammaLib distinguishes two types of events: event
atoms, which are individual events, and event bins, which are
collections of events with similar properties. Event atoms are
implemented by the abstract base class :doxy:`GEventAtom`, while event bins are
implemented by the abstract base class :doxy:`GEventBin`. Both classes derive
from the abstract :doxy:`GEvent` base class.

Each event type has it's own container class, which derives from the
abstract :doxy:`GEvents` base class. Event atoms are collected by the abstract
:doxy:`GEventList` base class, while event bins are collected by the abstract
:doxy:`GEventCube` base class. The :doxy:`GEventList` class contains an instance of the
abstract :doxy:`GRoi` base class.

The basic constitutents on an event are an energy, implemented by the :doxy:`GEnergy`
class, a time, implemented by the :doxy:`GTime` class and a direction, implemented by
the abstract :doxy:`GInstDir` base class. The latter is not necessarily a direction on
the sky, but can also be a detector number or any other information that encodes spatial
information for a given gamma-ray telescope.

The :doxy:`GEnergy` class implements energy values in a unit independent way, supporting
also conversion of energy values to specific units. The :doxy:`GTime` class accomplishes
the same for time, allowing also transformations to different time systems. Instances
of :doxy:`GTime` and :doxy:`GEnergy` can be collected using the :doxy:`GTimes` and
:doxy:`GEnergies` container classes. Time are also used to defined so called Good Time
Intervals, which are implemented by the :doxy:`GGti` class. Good Time Intervals are
specified in a given time reference, which is defined by the :doxy:`GTimeReference` class.

Instances of :doxy:`GEnergy`, :doxy:`GTime` and :doxy:`GSkyDir` (a direction on the sky)
compose the attributes of a photon, implemented by the :doxy:`GPhoton` class. Instances
of this class can be collected using the :doxy:`GPhotons` container class. Alternatively,
instead of a sky direction a spatial model can be used to characterise a gamma-ray source,
which is implemented using the :doxy:`GSource` class.

Support for pulsars is implemented by the :doxy:`GPulsar` class which is a collected of
ephemerides data, implemented by the :doxy:`GPulsarEphemeris` class.

Additional support classes of the module comprise the :doxy:`GPhases` class, which defines
phase intervals for pulsar or gamma-ray binary analysis, and the :doxy:`GEphemerides` class
which implements the JPL DE200 ephemerides that can be used for Barycentric corrections.

