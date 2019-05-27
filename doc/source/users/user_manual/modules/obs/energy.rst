.. _um_obs_energy:

Energies in GammaLib
====================

Energies are implemented by the :doxy:`GEnergy` class that transparently
handles energies with different units. The energy is stored internally in
MeV, but energies can also be provided in erg, keV, GeV, TeV, as well as
the equivalent wavelength in Angstrom.

The code threads below illustrate how an energy of 511 keV can be set using
the energy constructor or a dedicated method, and how that energy can be
retrieved in units of MeV.

**Python**

.. code-block:: python
   :linenos:

   import gammalib                                      # Make GammaLib available
   energy1 = gammalib.GEnergy(511.0,'keV')              # Using energy constructor
   energy2 = gammalib.GEnergy()
   energy2.keV(511.0)                                   # Using keV method
   E_MeV1 = energy2('MeV')                              # Get energy in MeV
   E_MeV2 = energy2.MeV()                               # Another variant to get energy in MeV

**C++**

.. code-block:: cpp
   :linenos:

   #include "GammaLib.hpp"                              // Make GammaLib available
   GEnergy energy1(511.0,"keV");                        // Using energy constructor
   GEnergy energy2;
   energy2.keV(511.0);                                  // Using keV method
   double E_MeV1 = energy2("MeV");                      // Get energy in MeV
   double E_MeV2 = energy2.MeV();                       // Another variant to get energy in MeV

Energies can be added and subtracted, or multiplied or divided by a floating
point value.
In addition, two energies can be divided, resulting in a floating point value.
Energies can also be compared using the usual operators.

Energies can be collected by the :doxy:`GEnergies` container class. The threads
below illustrate how energies can be appended or inserted into the
:doxy:`GEnergies` container class.

**Python**

.. code-block:: python
   :linenos:

   import gammalib                                      # Make GammaLib available
   energy1  = gammalib.GEnergy(511.0,'keV')             # Set first energy
   energy2  = gammalib.GEnergy(1024.0,'keV')            # Set second energy
   energies = gammalib.GEnergies()                      # Create energy container instance
   energies.append(energy1)                             # Append first energy
   energies.insert(0,energy2)                           # Insert second energy before first energy

**C++**

.. code-block:: cpp
   :linenos:

   #include "GammaLib.hpp"                              // Make GammaLib available
   GEnergy   energy1(511.0,"keV");                      // Set first energy
   GEnergy   energy2(1024.0,"keV");                     // Set second energy
   GEnergies energies;                                  // Create energy container instance
   energies.append(energy1);                            // Append first energy
   energies.insert(0,energy2);                          // Insert second energy before first energy

The energy container can be saved into a file using the :doxy:`GEnergies::save()`
method, and loaded from a file using either the load constructor or the
:doxy:`GEnergies::load()` method.

Another energy container class is :doxy:`GEbounds` which stores the energy
boundaries for an arbitrary number of energy bins.

