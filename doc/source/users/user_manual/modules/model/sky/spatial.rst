Spatial components
^^^^^^^^^^^^^^^^^^

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

An alternative XML format is supported for compatibility with the Fermi/LAT XML
format:

.. code-block:: xml

   <source name="Crab" type="PointSource">
     <spatialModel type="SkyDirFunction">
       <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
       <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
     </spatialModel>
     <spectrum type="...">
       ...
     </spectrum>
   </source>


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


Radial ring
===========

The ``RadialRing`` model specifies a uniform intensity distribution within
a circular ring. The circular ring is defined by the celestial coordinates
``RA`` and ``DEC`` of the ring centre, the ring **inner radius** defined by
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


Radial shell
============

The ``RadialShell`` model specifies a 3-dimensional shell projected on the
sky. The shell is defined by the celestial coordinates ``RA`` and ``DEC`` of
the shell centre, the **inner radius** of the shell defined by ``Radius`` and
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
