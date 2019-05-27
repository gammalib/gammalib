Under the hood
==============

.. _um_registry:

Registry classes
----------------

Registry classes statically hold instances of derived classes, and can be
used to allocate a derived class dependent on some text string.
The code thread below illustrates how the :doxy:`GModels::read()` method
extracts models from an XML file, allocates the requested model through the
:doxy:`GModelRegistry::alloc()` method, and appends the model to the model
container.

**C++**

.. code-block:: cpp
   :linenos:

   int n = lib->elements("source");                           // Get number of sources
   for (int i = 0; i < n; ++i) {                              // Loop over sources
       const GXmlElement* src = lib->element("source", i);    // Get pointer on source
       std::string type = src->attribute("type");             // Get source model type
       GModelRegistry registry;                               // Get model registry
       GModel*        ptr = registry.alloc(type);             // Allocate model of specific type
       if (ptr != NULL) {                                     // If model is valid:
           ptr->read(*src);                                   // - read model
           append(*ptr);                                      // - append it to container
           delete ptr;                                        // - free model
       }
   }

To access the registry information it is sufficient to create an instance of
the registry. To illustrate how the registry is filled the
``src/model/GModelSky.cpp`` file can be examined:

**C++**

.. code-block:: cpp
   :linenos:

   const GModelSky         g_pointsource_seed("PointSource");
   const GModelSky         g_extendedsource_seed("ExtendedSource");
   const GModelSky         g_diffusesource_seed("DiffuseSource");
   const GModelSky         g_compositesource_seed("CompositeSource");
   const GModelRegistry    g_pointsource_registry(&g_pointsource_seed);
   const GModelRegistry    g_extendedsource_registry(&g_extendedsource_seed);
   const GModelRegistry    g_diffusesource_registry(&g_diffusesource_seed);
   const GModelRegistry    g_compositesource_registry(&g_compositesource_seed);

Here four instances of :doxy:`GModelSky` are created as global variables which
are appended to the :doxy:`GModelRegistry` registry by creating an instance
of the registry which takes the address of the :doxy:`GModelSky` instance
as argument.

Registries are widely used throughout GammaLib for the handling of models and
its components, and the handling of the observations of the different supported
instruments.
