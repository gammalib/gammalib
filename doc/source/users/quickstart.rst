.. _quickstart:

Getting started
===============

.. _firststeps:

First steps with GammaLib
-------------------------

GammaLib comes with a Python interface, and as first step you should
verify that the Python interface works correctly. You do this by
typing:

.. code-block:: python
   :linenos:

   import gammalib
   gammalib.test()

The first command loads the GammaLib Python module, the second command 
executes the Python unit tests. If all tests are ok the GammaLib Python
module has been installed successfully. In case of problems, make sure 
that you've setup correctly the GammaLib environment. You basically need the
``PYTHONPATH`` variable set to GammaLib's Python site-package, and you also
have to make sure that GammaLib and any dependent library (cfitsio,
readline, ncurses) is in the library load path (eventually you may need
to adjust ``LD_LIBRARY_PATH``, or ``DYLD_LIBRARY_PATH`` if you are on Mac 
OS X). This is all done automatically if you set up the environment as 
described :ref:`here <setup_env>`.

Now try:

.. code-block:: python
   :linenos:

   models = gammalib.GModels()
   print(models)
   === GModels ===
   Number of models ..........: 0
   Number of parameters ......: 0

You just alloacted your first GammaLib object, which is an empty model
container.

Now let's append a model to this container. For this, type:

.. code-block:: python
   :linenos:

   pos = gammalib.GSkyDir()
   pos.radec_deg(83.6331, 22.0145)
   spatial = gammalib.GModelSpatialPointSource(pos)
   spectral = gammalib.GModelSpectralPlaw(1.0, -2.0, gammalib.GEnergy(100, 'MeV'))
   model = gammalib.GModelSky(spatial, spectral)
   models.append(model)
   print(models)
   === GModels ===
    Number of models ..........: 1
    Number of parameters ......: 6
   === GModelSky ===
    Name ......................:
    Instruments ...............: all
    Observation identifiers ...: all
    Model type ................: PointSource
    Model components ..........: "PointSource" * "PowerLaw" * "Constant"
    Number of parameters ......: 6
    Number of spatial par's ...: 2
     RA .......................: 83.6331 deg (fixed,scale=1)
     DEC ......................: 22.0145 deg (fixed,scale=1)
    Number of spectral par's ..: 3
     Prefactor ................: 1 +/- 0 [0,infty[ ph/cm2/s/MeV (free,scale=1,gradient)
     Index ....................: -2 +/- 0 [10,-10]  (free,scale=-2,gradient)
     PivotEnergy ..............: 100 MeV (fixed,scale=100,gradient)
    Number of temporal par's ..: 1
     Normalization ............: 1 (relative value) (fixed,scale=1,gradient)
    Number of scale par's .....: 0

With this sequence of commands you first defined a sky direction in
celestial coordinates using :doxy:`GSkyDir`. You may recognise that this
is the position of the Crab. You then used this position to define the
spatial component of a sky model using :doxy:`GModelSpatialPointSource`. As the name
suggests, the spatial component is a point source. Next, you defined the
spectral component using :doxy:`GModelSpactralPlaw`: a power law with
normalisation of 1 and a spectral index of -2. Then, you combined the spatial
and spectral components in a sky model using :doxy:`GModelSky`. And
finally you appended the sky model to the model container allocated
previously using the append method. When you then print the model
container you see that is contains now one model with 6 parameters.
Among them, you find the specified source position (parameters ``RA`` and
``DEC``), the power law normalisation (parameter ``Prefactor``) and the spectral
index (parameter ``Index``). In addition, the reference energy for the
normalisation has been set to 100 MeV (parameter ``PivotEnergy``)
and the temporal component has been set by default to a constant
(parameter ``Constant``).

Suppose you want to change one of these parameters, such as the
``PivotEnergy`` that you want to set to 1 TeV. This can be done using:

.. code-block:: python
   :linenos:

   models[0]['PivotEnergy'].value(1e6)
   print(models)
   ...
   PivotEnergy ..............: 1000000 MeV (fixed,scale=100,gradient)
   ...

As the units are MeV, we had to specify a value of 1e6 to set the
reference energy to 1 TeV. We did this by accessing the first model in
the container using ``models[0]`` (counting in GammaLib always starts from
0). Then we addressed the ``PivotEnergy`` parameter by specifying
``['PivotEnergy']``. And finally we called the value method that sets the
value of a particular parameter.

After all this hard work, you may save your model into a XML file by
typing:

.. code-block:: python
   :linenos:

   models.save('test.xml')

and you can load it from an XML file in memory using:

.. code-block:: python
   :linenos:

   new_models = gammalib.GModels('test.xml')
   print(new_models)

The last print command is to convince yourself that the models have been
loaded properly.

Much more is still to come. Please be a little bit patient, we're working
on it. In the meantime you may check the `Doxygen
documentation <../doxygen/index.html>`_ to see what classes and methods are
available.

Getting Help
------------

Any questions, bug reports, or suggested enhancements related to
GammaLib should be submitted via the
`issue tracker <https://cta-redmine.irap.omp.eu/projects/gammalib/issues/new>`_
or the
`mailing list <mailto:ctools@irap.omp.eu>`_.
