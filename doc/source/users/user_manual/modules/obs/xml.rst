Describing observations using XML
=================================

Observations can be described in GammaLib using an ASCII file in XML
format. The general format of this file is

.. code-block:: xml

    <observation_list title="observation library">
      <observation name="..." id="..." instrument="...">
        ...
      </observation>
      <observation name="..." id="..." instrument="...">
        ...
      </observation>
      ...
    </observation_list>

where each ``<observation>`` tag describes one observation. Each observation
has a ``name`` attribute, an ``id`` (identifier) attribute and an
``instrument`` attribute. The latter decides which instrument specific
class will be allocated upon reading the XML file. For a given instrument,
observation identifiers must be unique.

The specific format of the XML file for a given instrument is defined by the
relevant instrument specific :doxy:`GObservation` class. For example, a CTA
observation implemented by the :doxy:`GCTAObservation` class is described by

.. code-block:: xml

    <observation name="..." id="..." instrument="...">
      <parameter name="EventList"           file="..."/>
      <parameter name="EffectiveArea"       file="..."/>
      <parameter name="PointSpreadFunction" file="..."/>
      <parameter name="EnergyDispersion"    file="..."/>
    </observation>

for an unbinned observation and by

.. code-block:: xml

    <observation name="..." id="..." instrument="...">
      <parameter name="CountsCube"          file="..."/>
      <parameter name="EffectiveArea"       file="..."/>
      <parameter name="PointSpreadFunction" file="..."/>
      <parameter name="EnergyDispersion"    file="..."/>
    </observation>

for a binned observation. Here, ``EventList`` specifies a FITS file containing
an event list and ``CountsCube`` specifies a FITS file containing a counts 
cube.
The other tags specify the components of the instrumental response function
and are optional.
Similar definitions exist for the other instruments.

The code threads below show different variants of how observations can be
loaded from or saved into an XML file.

**Python**

.. code-block:: python

   >>> import gammalib                                      # Make GammaLib available
   >>> obs = gammalib.GObservations('my_observations.xml')  # Construct observations from XML file
   >>> obs.save('my_copied_observations.xml')               # Save observations into XML file

.. code-block:: python

   >>> import gammalib                                      # Make GammaLib available
   >>> obs = gammalib.GObservations()                       # Construct empty observations
   >>> obs.load('my_observations.xml')                      # Load observations from XML file
   >>> obs.save('my_copied_observations.xml')               # Save observations into XML file

.. code-block:: python

   >>> import gammalib                                      # Make GammaLib available
   >>> obs    = gammalib.GObservations()                    # Construct empty observations
   >>> xmlin  = gammalib.GXml('my_observations.xml')        # Construct GXml instance from XML file
   >>> xmlout = gammalib.GXml()                             # Construct empty GXml instance
   >>> obs.read(xmlin)                                      # Read observations from GXml instance
   >>> obs.write(xmlout)                                    # Write observations into GXml instance
   >>> xmlout.save('my_copied_observations.xml')            # Save GXml instance in XML file

**C++**

.. code-block:: cpp

   #include "GammaLib.hpp"                                  // Make GammaLib available
   GObservations obs("my_observations.xml");                // Construct observations from XML file
   obs.save("my_copied_observations.xml");                  // Save observations into XML file

.. code-block:: cpp

   #include "GammaLib.hpp"                                  // Make GammaLib available
   GObservations obs;                                       // Construct empty observations
   obs.load("my_observations.xml");                         // Load observations from XML file
   obs.save("my_copied_observations.xml");                  // Save observations into XML file

.. code-block:: cpp

   #include "GammaLib.hpp"                                  // Make GammaLib available
   GObservations obs;                                       // Construct empty observations
   GXml          xmlin("my_observations.xml");              // Construct GXml instance from XML file
   GXml          xmlout;                                    // Construct empty GXml instance
   obs.read(xmlin);                                         // Read observations from GXml instance
   obs.write(xmlout);                                       // Write observations into GXml instance
   xmlout.save("my_copied_observations.xml");               // Save GXml instance in XML file
