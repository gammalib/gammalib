XML model definition format
===========================

Sky models can be defined using an XML file format that is inspired by the format
used for the Large Area Telescope (LAT) aboard NASA's Fermi satellite. To
illustrate the format, the definition of a point source at the position of
the Crab nebula with a power spectral shape is shown below. The tags
``<spatialModel>``, ``<spectrum>`` and ``<temporalModel>`` describe the
spatial, spectral and temporal model components defined above. The
``<temporalModel>`` tag can be omitted, which will automatically allocate a
constant temporal component for the model.

Each source model needs to have a unique name specified by the ``name``
attributes. The ``type`` attribute specifies the type of the model and
needs to be one of ``PointSource``, ``ExtendedSource``, ``DiffuseSource``
and ``CompositeSource``. Optional attributes include:

* ``id``, specifying an identifier string for a model that can be used to tie
  the source model to one or several observations;,
* ``ts``, specifying the Test Statistic value of the source model component, and
* ``tscalc``, specifying whether the Test Statistic value should be computed
  in a maximum-likelihood computation.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      <spatialModel type="PointSource">
        <parameter name="RA"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
        <parameter name="DEC" scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
      </spatialModel>
      <spectrum type="PowerLaw">
        <parameter name="Prefactor"   scale="1e-16" value="5.7"  min="1e-07" max="1000.0" free="1"/>
        <parameter name="Index"       scale="-1"    value="2.48" min="0.0"   max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"   value="0.3"  min="0.01"  max="1000.0" free="0"/>
      </spectrum>
      <temporalModel type="Constant">
        <parameter name="Normalization" scale="1.0" value="1.0" min="0.1" max="10.0" free="0"/>
      </temporalModel>
    </source>
  </source_library>

To accomodate for possible calibration biases for different instruments,
instrument specific scaling factors can be specified as part of a model.
In the example below, the model will be multiplied by a factor of 1 for ``LAT``
and a factor of 0.5 for ``CTA``.

.. code-block:: xml

  <?xml version="1.0" standalone="no"?>
  <source_library title="source library">
    <source name="Crab" type="PointSource">
      ...
      <scaling>
        <instrument name="LAT" scale="1.0" min="0.1" max="10.0" value="1.0" free="0"/>
        <instrument name="CTA" scale="1.0" min="0.1" max="10.0" value="0.5" free="0"/>
      </scaling>
    </source>
  </source_library>
