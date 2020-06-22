Spatial components
^^^^^^^^^^^^^^^^^^

Point source
============

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


Radial ring
===========

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


Radial Gaussian
===============

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
