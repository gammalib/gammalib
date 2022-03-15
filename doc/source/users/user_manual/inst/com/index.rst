.. um_com:

COMPTEL interface
-----------------

Overview
~~~~~~~~

The COMPTEL module provides an instrument interface for the COMPTEL 
telescope that has been operated on the Compton Gamma-Ray Observatory 
(CGRO) from 1991 to 2000.

The following figure presents an overview over the classes of the COMPTEL
module and their relations.

.. _fig_uml_com:

.. figure:: uml_com.png
   :width: 90%
   :align: center

   *COMPTEL module*

The central object is a COMPTEL observation, implemented by the
:doxy:`GCOMObservation` class. The class either holds a COMPTEL event cube,
implemented by the :doxy:`GCOMEventCube` class, or an event list, implemented
by the :doxy:`GCOMEventList` class. A single bin of the event cube is
implemented by the :doxy:`GCOMEventBin` class, a single event by the
:doxy:`GCOMEventAtom` class.

The :doxy:`GCOMEventCube` class uses the general :doxy:`GCOMDri` class to
store the data. :doxy:`GCOMDri` implements a general three-dimensional
data-space of COMPTEL, spanned by the scatter direction :math:`\chi` and
:math:`\Psi` and the scatter angle :math:`\bar{\varphi}`.
This three-dimensional coordinate is implemented by the :doxy:`GCOMInstDir`
class.
The :doxy:`GCOMDri` class further stores the selection sets, implemented by the
:doxy:`GCOMSelection` class.

The :doxy:`GCOMObservation` class further holds COMPTEL Good Time Intervals,
implemented by the :doxy:`GCOMTim` class, and Orbit Aspect Data, implemented
by the :doxy:`GCOMOads` container class that holds records of :doxy:`GCOMOad`.
Optionally, the data may also hold information needed for barycentre correction
of event arrival times, implemented by the :doxy:`GCOMBvcs` container class that
holds records of :doxy:`GCOMBvc`. If this information is not provided, the
barycentre correction information is computed on the fly, using JPL DE200
ephemerides that are provided by the :doxy:`GEphemerides` class.

The module also provides a :doxy:`GCOMIaq` class that enables computation of
instrument response functions. The computation is based on the D1 and D2
detector response functions, that are accessed via the :doxy:`GCOMD1Response`
and :doxy:`GCOMD2Response` classes. Instrument characteristics that are also
required for the computation are access via the :doxy:`GCOMInstChars` class.

So far two COMPTEL background model exists, that performs fitting of the
:math:`\bar{\varphi}` layers of a DRI model cube. While the
:doxy:`GCOMModelDRBPhibarBins` class fits each :math:`\bar{\varphi}` layer
with a separate scaling factor, the :doxy:`GCOMModelDRBPhibarNodes` class
specifies a scaling factor for a number of :math:`\bar{\varphi}` values and
interpolates the results linearly for :math:`\bar{\varphi}` values that are
not represented by the nodes. In addition, the :doxy:`GCOMModelDRM` class implements
the fitting of a data space model with a single scaling factor, and may be used
either for source or background model fitting.

Finally, the :doxy:`GCOMStatus` class holds instrument status information.
