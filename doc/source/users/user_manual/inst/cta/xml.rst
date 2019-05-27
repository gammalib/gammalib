.. _um_cta_xml:

Describing CTA observations using XML
=====================================

CTA observations can be described in GammaLib using an ASCII file in XML
format (see :ref:`um_obs_xml`). The CTA specific section of this file has
the format:

.. code-block:: xml

   <observation name="..." id="..." instrument="...">
       <parameter name="EventList"           file="..."/>
       <parameter name="EffectiveArea"       file="..."/>
       <parameter name="PointSpreadFunction" file="..."/>
       <parameter name="EnergyDispersion"    file="..."/>
       <parameter name="Background"          file="..."/>
   </observation>

for an unbinned observation and

.. code-block:: xml

   <observation name="..." id="..." instrument="...">
       <parameter name="CountsCube"          file="..."/>
       <parameter name="EffectiveArea"       file="..."/>
       <parameter name="PointSpreadFunction" file="..."/>
       <parameter name="EnergyDispersion"    file="..."/>
       <parameter name="Background"          file="..."/>
   </observation>
 
for a binned observation.
Each parameter within the ``<observation>`` tag specifies the filename
for a specific file that is needed for the analysis.
The difference between an unbinned and a binned observation is that the 
first uses an event list (requiring a parameter with name ``EventList``)
while the latter uses an event cube  (requiring a parameter with name
``CountsCube``).
The ``EffectiveArea``, ``PointSpreadFunction``, ``EnergyDispersion`` and
``Background`` parameters provide the filenames of the instrument
specific response function components for each observation 
(see :ref:`um_cta_response`).
Alternatively, the calibration database and response name can be
specified, e.g.

.. code-block:: xml

   <observation name="..." id="..." instrument="...">
       <parameter name="EventList" file="..."/>
       <parameter name="Calibration" database="..." response="..."/>
   </observation>

A variant of the binned analysis is the stacked analysis.
While in a binned analysis an event cube is generated for each 
observation, event cubes will be summed in a stacked analysis.
Summing the events requires computation of the total exposure,
the average point spread function and the average background
rate (energy dispersion is not yet handled for stacked analysis).
Information for a stacked observation is provided in the
following format:

.. code-block:: xml

   <observation name="..." id="..." instrument="...">
       <parameter name="CountsCube"   file="..."/>
       <parameter name="ExposureCube" file="..."/>
       <parameter name="PsfCube"      file="..."/>
       <parameter name="BkgCube"      file="..."/>
   </observation>

The stacked analysis uses also an event cube but now requires the
pre-computed response cubes.
The ``ExposureCube``, ``PsfCube``, and ``BkgCube`` parameters
provide the filenames of the required components.

The ``instrument`` attribute of CTA observations can be one of
``CTA``, ``HESS``, ``MAGIC`` or ``VERITAS``. This allows mixing of
observations from difference IACTs within a single analysis.
Note that no code that is specific to any of these four instruments is
implemented in GammaLib, but the ``instrument`` attribute is used to
tie models to instruments, allowing thus to provide specific background
models for each of the instruments in a combined analysis. 

The ``id`` attribute specifies an identifier for the observation that 
needs to be unique for a given instrument.
The identifier allows to connect an observation to a specific model
component.

Optionally, an observation can also have user defined minimum and maximum 
energy boundaries.
These boundaries are specified as optional attributes in the XML file, e.g. 

.. code-block:: xml

   <observation name="..." id="..." instrument="..." emin="0.1" emax="100.0">
       ...
   </observation>

.. note::
   Energy boundaries in the observation XML file are specified in units
   of TeV.

In case that neither an event list nor an event cube is available, a
CTA observation can be defined by specifying the pointing direction,
the Good Time Intervals, the region of interest, and optionally the
deadtime correction factor.
Here an example of the expected XML format:

.. code-block:: xml

   <observation name="GPS" id="000001" instrument="CTA">
       <parameter name="Pointing" ra="186.721" dec="-61.4328" />
       <parameter name="GoodTimeIntervals" tmin="0" tmax="35100" />
       <parameter name="TimeReference" mjdrefi="51544" mjdreff="0.5" timeunit="s" timesys="TT" timeref="LOCAL" />
       <parameter name="RegionOfInterest" ra="186.721" dec="-61.4328" rad="5" />
       <parameter name="Deadtime" deadc="0.95" />
       <parameter name="Calibration" database="prod2" response="South_50h" />
   </observation>

.. note::
   The time reference for the Good Time Intervals is specified 
   using the ``TimeReference`` parameter.
