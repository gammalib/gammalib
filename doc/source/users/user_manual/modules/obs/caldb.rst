.. _um_obs_caldb:

Using a calibration database
============================

The class :doxy:`GCaldb` implements the interface to a calibration database
that is compliant with the HEASARC format.
To use a calibration database, make sure that the environment variable
``CALDB`` is set to a valid calibration database.
The database is the open by specifying a mission and an instrument name.
The examples below show two variants of opening the ``prod2`` database for
``CTA``.

**Python**

.. code-block:: python
   :linenos:

   import gammalib                                      # Make GammaLib available
   caldb1 = gammalib.GCaldb('cta','prod2')              # Using opening constructor
   caldb2 = gammalib.GCaldb()
   caldb2.open('cta','prod2')                           # Using open method

**C++**

.. code-block:: cpp
   :linenos:

   #include "GammaLib.hpp"                              // Make GammaLib available
   GCaldb caldb1("cta","prod2");                        // Using opening constructor
   GCaldb caldb2;
   caldb2.open("cta","prod2");                          // Using open method

Once opened, the database can be used to locate for example the filename
of a response component. The examples below show how the effective area
filename for the response ``South_5h`` is accessed via the calibration
database.

**Python**

.. code-block:: python
   :linenos:

   filename = caldb1.filename('','','EFF_AREA','','','NAME(South_5h)')

**C++**

.. code-block:: cpp
   :linenos:

   GFilename filename = caldb1.filename("","","EFF_AREA","","","NAME(South_5h)");
