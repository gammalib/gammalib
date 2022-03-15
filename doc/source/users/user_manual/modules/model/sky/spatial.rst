Spatial components
^^^^^^^^^^^^^^^^^^

.. note::
   In the following model descriptions, celestial coordinates ``RA`` and ``DEC``
   may be replaced by Galactic coordinates ``GLON`` and ``GLAT``.


Point source
============

The ``PointSource`` model specifies a source that has no spatial extension.
It is defined by its celestial coordinates ``RA`` and ``DEC`` given in
units of degrees.

.. code-block:: xml

   <source name="Crab" type="PointSource">
     <spatialModel type="PointSource">
       <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

where

* ``RA`` is the Right Ascension (degrees)
* ``DEC`` is the Declination (degrees)

.. note::
   For compatibility with the Fermi/LAT ScienceTools the model type
   ``PointSource`` can be replaced by ``SkyDirFunction``.


Radial disk
===========

The ``RadialDisk`` model specifies a uniform circular intensity distribution,
defined by the celestial coordinates ``RA`` and ``DEC`` of the disk centre
and the disk ``Radius``. All parameters are given in units of degrees.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="RadialDisk">
       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
       <parameter name="Radius" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

where

* ``RA`` is the Right Ascension of the disk centre (degrees)
* ``DEC`` is the Declination of the disk centre (degrees)
* ``Radius`` is the disk radius (degrees)


Radial ring
===========

The ``RadialRing`` model specifies a uniform intensity distribution within
a circular ring. The circular ring is defined by the celestial coordinates
``RA`` and ``DEC`` of the ring centre, the ring inner radius defined by
``Radius`` and the ring width, defined by ``Width``. Specifically, the
ring outer radius is given by ``Radius+Width``. All parameters are given
in units of degrees.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="RadialRing">
       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
       <parameter name="Radius" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
       <parameter name="Width"  scale="1.0" value="0.15"    min="0.01" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

where

* ``RA`` is the Right Ascension of the ring centre (degrees)
* ``DEC`` is the Declination of the ring centre (degrees)
* ``Radius`` is the inner ring radius (degrees)
* ``Width`` is the ring width radius (degrees)


Radial Gaussian
===============

The ``RadialGaussian`` model specifies a spherical Gaussian intensity
distribution, defined by the celestial coordinates ``RA`` and ``DEC`` of the
Gaussian centre and the Gaussian ``Sigma`` parameter. All parameters are given
in units of degrees.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="RadialGaussian">
       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
       <parameter name="Sigma" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

and implements

.. math::
   M_{\rm spatial}(\theta) = \frac{1}{2 \pi \sigma^2} \exp
                  \left(-\frac{1}{2}\frac{\theta^2}{\sigma^2} \right),

where

* ``RA`` is the Right Ascension of the Gaussian centre (degrees)
* ``DEC`` is the Declination of the Gaussian centre (degrees)
* :math:`\sigma` = ``Sigma`` (degrees)


Radial general Gaussian
=======================

The ``RadialGeneralGaussian`` model specifies a generalised Gaussian intensity distribution,
defined by the celestial coordinates ``RA`` and ``DEC`` of the generalised Gaussian centre,
a radius ``Radius``and a radial index parameter ``R_Index``.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="RadialGeneralGaussian">
       <parameter name="RA"      scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC"     scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
       <parameter name="Radius"  scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
       <parameter name="R_Index" scale="1.0" value="0.5"     min="0.01" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

and implements

.. math::
   M_{\rm spatial}(\theta) = \frac{1}{2 \pi r^2 \eta \Gamma(2\eta)} \exp
                  \left[- \left(\frac{\theta^2}{r^2}\right)^\frac{1}{\eta} \right],

where

* ``RA`` is the Right Ascension of the Gaussian centre (degrees)
* ``DEC`` is the Declination of the Gaussian centre (degrees)
* :math:`r` = ``Radius`` (degrees)
* :math:`\eta` = ``R_Index``

The model normalisation is correct in the small angle approximation and for
:math:`\eta` of the order of unity or smaller.


Radial shell
============

The ``RadialShell`` model specifies a 3-dimensional shell projected on the
sky. The shell is defined by the celestial coordinates ``RA`` and ``DEC`` of
the shell centre, the inner radius of the shell defined by ``Radius`` and
the width of the shell, defined by ``Width``. Specifically, the outer radius
of the shell is given by ``Radius+Width``. All parameters are given in units
of degrees.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="RadialShell">
       <parameter name="RA"     scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC"    scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
       <parameter name="Radius" scale="1.0" value="0.30"    min="0.01" max="10"  free="1"/>
       <parameter name="Width"  scale="1.0" value="0.10"    min="0.01" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

and implements

.. math::
   M_{\rm spatial}(\theta) =  n_0 \left \{
   \begin{array}{l l}
      \displaystyle
      \sqrt{ \theta_{\rm out}^2 - \theta^2 } - \sqrt{ \theta_{\rm in}^2 - \theta^2 }
      & \mbox{if $\theta \le \theta_{\rm in}$} \\
      \\
     \displaystyle
      \sqrt{ \theta_{\rm out}^2 - \theta^2 }
      & \mbox{if $\theta_{\rm in} < \theta \le \theta_{\rm out}$} \\
      \\
     \displaystyle
     0 & \mbox{if $\theta > \theta_{\rm out}$}
   \end{array}
   \right .

where

* ``RA`` is the Right Ascension of the shell centre (degrees)
* ``DEC`` is the Declination of the shell centre (degrees)
* :math:`\theta_{\rm out}` = ``Radius`` + ``Width`` (degrees)
* :math:`\theta_{\rm in}` = ``Radius`` (degrees)


Radial profiles
===============

Radial profiles are defined by a arbitrary function of the radial distance from a central
position. The following radial profiles exist:

Burkert Dark matter profile
---------------------------

.. code-block:: xml

  <source name="Crab" type="ExtendedSource">
    <spatialModel type="DMBurkertProfile">
      <parameter name="RA"           scale="1.0" value="83.6331" min="-360"    max="360"   free="1"/>
      <parameter name="DEC"          scale="1.0" value="22.0145" min="-90"     max="90"    free="1"/>
      <parameter name="ScaleRadius"  scale="1.0" value="21.5"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="ScaleDensity" scale="1.0" value="0.2"     min="0.0001"  max="1000"  free="0"/>
      <parameter name="HaloDistance" scale="1.0" value="7.94"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="ThetaMin"     scale="1.0" value="1.0e-6"  min="1.0e-10" max="1000"  free="0"/>
      <parameter name="ThetaMax"     scale="1.0" value="180.0"   min="0.0001"  max="1000"  free="0"/>
      <parameter name="CoreRadius"   scale="1.0" value="0.5"     min="0.0001"  max="1000"  free="0"/>
    </spatialModel>
    <spectrum type="...">
      ...
    </spectrum>
  </source>

Einasto Dark matter profile
---------------------------

.. code-block:: xml

  <source name="Crab" type="ExtendedSource">
    <spatialModel type="DMEinastoProfile">
      <parameter name="RA"           scale="1.0" value="83.6331" min="-360"    max="360"   free="1"/>
      <parameter name="DEC"          scale="1.0" value="22.0145" min="-90"     max="90"    free="1"/>
      <parameter name="ScaleRadius"  scale="1.0" value="21.5"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="ScaleDensity" scale="1.0" value="0.2"     min="0.0001"  max="1000"  free="0"/>
      <parameter name="HaloDistance" scale="1.0" value="7.94"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="Alpha"        scale="1.0" value="0.17"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="ThetaMin"     scale="1.0" value="1.0e-6"  min="1.0e-10" max="1000"  free="0"/>
      <parameter name="ThetaMax"     scale="1.0" value="180.0"   min="0.0001"  max="1000"  free="0"/>
      <parameter name="CoreRadius"   scale="1.0" value="0.5"     min="0.0001"  max="1000"  free="0"/>
    </spatialModel>
    <spectrum type="...">
      ...
    </spectrum>
  </source>

Zhao Dark matter profile
------------------------

.. code-block:: xml

  <source name="Crab" type="ExtendedSource">
    <spatialModel type="DMZhaoProfile">
      <parameter name="RA"           scale="1.0" value="83.6331" min="-360"    max="360"   free="1"/>
      <parameter name="DEC"          scale="1.0" value="22.0145" min="-90"     max="90"    free="1"/>
      <parameter name="ScaleRadius"  scale="1.0" value="21.5"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="ScaleDensity" scale="1.0" value="0.2"     min="0.0001"  max="1000"  free="0"/>
      <parameter name="HaloDistance" scale="1.0" value="7.94"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="Alpha"        scale="1.0" value="0.17"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="Beta"         scale="1.0" value="3.00"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="Gamma"        scale="1.0" value="1.00"    min="0.0001"  max="1000"  free="0"/>
      <parameter name="ThetaMin"     scale="1.0" value="1.0e-6"  min="1.0e-10" max="1000"  free="0"/>
      <parameter name="ThetaMax"     scale="1.0" value="180.0"   min="0.0001"  max="1000"  free="0"/>
      <parameter name="CoreRadius"   scale="1.0" value="0.5"     min="0.0001"  max="1000"  free="0"/>
    </spatialModel>
    <spectrum type="...">
      ...
    </spectrum>
  </source>

Gaussian profile
----------------

This profile is equivalent to ``RadialGaussian``.

.. code-block:: xml

  <source name="Crab" type="ExtendedSource">
    <spatialModel type="GaussianProfile">
      <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
      <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
      <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
    </spatialModel>
    <spectrum type="...">
      ...
    </spectrum>
  </source>


Elliptical disk
===============

The ``EllipticalDisk`` model specifies a uniform elliptical intensity
distribution, defined by the celestial coordinates ``RA`` and ``DEC`` of the
centre of the ellipse, the minor and major radii ``MinorRadius`` and
``MajorRadius`` of the ellipse, and the position angle ``PA`` that is
counted counter-clockwise from celestial North. All parameters are given in
units of degrees.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="EllipticalDisk">
       <parameter name="RA"          scale="1.0" value="83.6331" min="-360"  max="360" free="1"/>
       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"   max="90"  free="1"/>
       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

where

* ``RA`` is the Right Ascension (degrees)
* ``DEC`` is the Declination (degrees)
* ``PA`` is the position angle, counted counterclockwise from North (degrees)
* ``MinorRadius`` is the minor radius of the ellipse (degrees)
* ``MajorRadius`` is the major radius of the ellipse (degrees)


Elliptical Gaussian
===================

The ``EllipticalGaussian`` model specifies an elliptical Gaussian intensity
distribution, defined by the celestial coordinates ``RA`` and ``DEC`` of the
centre of the ellipse, the minor and major sigma parameter ``MinorRadius`` and
``MajorRadius`` of the ellipse, and the position angle ``PA`` that is
counted counter-clockwise from celestial North. All parameters are given in
units of degrees.

.. code-block:: xml

   <source name="Crab" type="ExtendedSource">
     <spatialModel type="EllipticalGaussian">
       <parameter name="RA"          scale="1.0" value="83.6331" min="-360"  max="360" free="1"/>
       <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"   max="90"  free="1"/>
       <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
       <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
       <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

and implements

.. math::
   M_{\rm spatial}(\theta, \phi) = \exp \left( -\frac{\theta^2}{2 r_\mathrm{eff}^2} \right),

with

.. math::
   r_\mathrm{eff} = \frac{ab} {\sqrt{\left( a \sin (\phi - \phi_0) \right)^2 +
                    \sqrt{\left( b \cos (\phi - \phi_0) \right)^2}}}

where

* ``RA`` is the Right Ascension (degrees)
* ``DEC`` is the Declination (degrees)
* ``PA`` is the position angle, counted counterclockwise from North (degrees)
* :math:`a` = ``MinorRadius`` (degrees)
* :math:`b` = ``MajorRadius`` (degrees)
* :math:`\phi_0` is the position angle of the ellipse, counted counterclockwise
  from North
* :math:`\phi` is the azimuth angle with respect to North.


EllipticalGeneralGaussian
=========================

The ``EllipticalGeneralGaussian`` model describes a Gaussian intensity distribution

.. code-block:: xml

  <source name="Crab" type="ExtendedSource">
    <spatialModel type="EllipticalGeneralGaussian">
      <parameter name="RA"          scale="1.0" value="83.6331" min="-360"  max="360" free="1"/>
      <parameter name="DEC"         scale="1.0" value="22.0145" min="-90"   max="90"  free="1"/>
      <parameter name="PA"          scale="1.0" value="45.0"    min="-360"  max="360" free="1"/>
      <parameter name="MinorRadius" scale="1.0" value="0.5"     min="0.001" max="10"  free="1"/>
      <parameter name="MajorRadius" scale="1.0" value="2.0"     min="0.001" max="10"  free="1"/>
      <parameter name="R_Index"     scale="1.0" value="0.5"     min="0.01"  max="10"  free="1"/>
    </spatialModel>
    <spectrum type="...">
      ...
    </spectrum>
  </source>

and implements

.. math::
   M_{\rm spatial}(\theta, \phi) = \frac{1}{2 \pi r^2 \eta
   \Gamma(2\eta)} \exp \left[ -\left(\frac{\theta^2}{2 r_\mathrm{eff}^2}\right)^\frac{1}{\eta} \right],

with

.. math::
   r_\mathrm{eff} = \frac{ab} {\sqrt{\left( a \sin (\phi - \phi_0) \right)^2 +
                    \sqrt{\left( b \cos (\phi - \phi_0) \right)^2}}}

where

* ``RA`` is the Right Ascension (degrees)
* ``DEC`` is the Declination (degrees)
* ``PA`` is the position angle, counted counterclockwise from North (degrees)
* :math:`a` = ``MinorRadius`` (degrees)
* :math:`b` = ``MajorRadius`` (degrees)
* :math:`\phi_0` is the position angle of the ellipse, counted counterclockwise
  from North
* :math:`\phi` is the azimuth angle with respect to North
* :math:`\eta` = ``R_Index``

The model normalisation is correct in the small angle approximation and for
:math:`\eta` of the order of unity or smaller.


Isotropic source
================

The ``DiffuseIsotropic`` model specifies an isotropic intensity distribution.
The only parameter of the model is a normalisation factor, specified by the
parameter ``Value``.

.. code-block:: xml

   <source name="Crab" type="DiffuseSource">
     <spatialModel type="DiffuseIsotropic">
       <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

An alternative XML format is supported for compatibility with the Fermi/LAT XML
format:

.. code-block:: xml

   <source name="Crab" type="DiffuseSource">
     <spatialModel type="ConstantValue">
       <parameter name="Value" scale="1" value="1" min="1"  max="1" free="0"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>


Diffuse map
===========

The ``DiffuseMap`` model specifies an intensity distribution that is
represented by a FITS image. The name of the FITS file is specified using
the ``file`` attribute of the ``spatialModel`` tag. If there are several
image in the FITS file, the first image will be extracted for the diffuse
map. Alternatively, the name of the relevant image extension or the extension
number can be specified in square brackets to select a specific image from
the FITS file.

The only parameter of the model is a normalisation factor, specified by the
parameter ``Normalization``.

.. code-block:: xml

   <source name="Crab" type="DiffuseSource">
     <spatialModel type="DiffuseMap" file="map.fits">
       <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

An alternative XML format is supported for compatibility with the Fermi/LAT XML
format:

.. code-block:: xml

   <source name="Crab" type="DiffuseSource">
     <spatialModel type="SpatialMap" file="map.fits">
       <parameter name="Prefactor" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
     </spatialModel>
     <spectrum type="...">
      ...
     </spectrum>
   </source>


Diffuse map cube
================

The ``DiffuseMapCube`` model specifies an energy-dependent intensity
distribution that is represented by a FITS file. The name of the FITS file is
specified using the ``file`` attribute of the ``spatialModel`` tag. The model
expects a 3-dimensional FITS image plus an extension with the name
``ENERGIES`` that specifies the energy for every layer of the FITS image.
The number of energies must correspond to the length of the 3rd image axis.

The only parameter of the model is a normalisation factor, specified by the
parameter ``Normalization``.

.. code-block:: xml

   <source name="Crab" type="DiffuseSource">
     <spatialModel type="DiffuseMapCube" file="map_cube.fits">
       <parameter name="Normalization" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>

An alternative XML format is supported for compatibility with the Fermi/LAT XML
format:

.. code-block:: xml

   <source name="Crab" type="DiffuseSource">
     <spatialModel type="MapCubeFunction" file="map_cube.fits">
       <parameter name="Value" scale="1" value="1" min="0.001" max="1000.0" free="0"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>


Composite model
===============

Spatial model components can be combined into a single model using the
:doxy:`GModelSpatialComposite class`. The class computes

.. math::
   M_{\rm spatial}(p|E,t) = \frac{1}{N} \sum_{i=0}^{N-1} M_{\rm spatial}^{(i)}(p|E,t)

where :math:`M_{\rm spatial}^{(i)}(p|E,t)` is any spatial model component
(including another composite model), and :math:`N` is the number of
model components that are combined.

An example of an XML file for a composite spatial model is shown below. In
this example, a point source is added to a radial Gaussian source to form
a composite spatial model. All spatial parameters of the composite model are
fitted.

.. code-block:: xml

   <source name="Crab" type="CompositeSource">
     <spatialModel type="Composite">
       <spatialModel type="PointSource" component="PointSource">
         <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
         <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
       </spatialModel>
       <spatialModel type="RadialGaussian">
         <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
         <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
         <parameter name="Sigma" scale="1.0" value="0.20"    min="0.01" max="10"  free="1"/>
       </spatialModel>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>
