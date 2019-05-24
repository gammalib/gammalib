Model factorisation
===================

Sky models describe the spatial, spectral and temporal properties of a 
gamma-ray source using the following factorisation:

.. math::
   M(p,E,t) = M_{\rm spatial}(p | E,t) \times
              M_{\rm spectral}(E | t) \times
              M_{\rm temporal}(t)

(:math:`M(p,E,t)` is given in units of
:math:`photons \,\, cm^{-2} s^{-1} MeV^{-1} sr^{-1}`).
The spatial model component :math:`M_{\rm spatial}(p|E,t)`
is defined by the abstract :doxy:`GModelSpatial` class, the spectral
model component :math:`M_{\rm spectral}(E|t)` is defined by the
abstract :doxy:`GModelSpectral` class and the temporal component
:math:`M_{\rm temporal}(t)` is defined by the abstract
:doxy:`GModelTemporal` class.

The spatial model component describes the energy and time dependent
morphology of the source.
It satisfies

.. math::
   \int_{\Omega} M_{\rm spatial}(p|E,t) \, d\Omega = 1

for all :math:`E` and :math:`t`, hence the spatial component does not
impact the spatially integrated spectral and temporal properties of the
source (the integration is done here over the spatial parameters
:math:`p` in a spherical coordinate system).
The units of the spatial model component are
:math:`[M_{\rm spatial}] = {\rm sr}^{-1}`.

The spectral model component describes the spatially integrated time
dependent spectral distribution of the source.
It satisfies

.. math::
   \int_{E} M_{\rm spectral}(E | t) \, dE = \Phi

for all :math:`t`, where :math:`\Phi` is the spatially and spectrally
integrated total source flux. The spectral component does not impact
the temporal properties of the integrated flux :math:`\Phi`.
The units of the spectral model component
are :math:`[M_{\rm spectral}] = {\rm cm}^{-2} {\rm s}^{-1} {\rm MeV}^{-1}`.

The temporal model component describes the relative variation of the
source flux with respect to the mean value given by the spectral model
component.
The temporal model component is unit less.
