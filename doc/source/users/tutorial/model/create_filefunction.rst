How to create a spectral file function?
=======================================

The following example illustrates how you can create a spectral file function
composed of four nodes and how you can save this file function into the file
``my_file_function.txt``.

**C++**

.. code-block:: cpp
   :linenos:

   GModelSpectralFunc spectrum;
   spectrum.append(GEnergy(1.0, "MeV"), 9.0e-5);
   spectrum.append(GEnergy(2.0, "MeV"), 7.0e-5);
   spectrum.append(GEnergy(4.0, "MeV"), 5.0e-5);
   spectrum.append(GEnergy(8.0, "MeV"), 3.0e-5);
   spectrum.save("my_file_function.txt", true);

**Python**

.. code-block:: python
   :linenos:

   spectrum = gammalib.GModelSpectralFunc()
   spectrum.append(gammalib.GEnergy(1.0, 'MeV'), 9.0e-5)
   spectrum.append(gammalib.GEnergy(2.0, 'MeV'), 7.0e-5)
   spectrum.append(gammalib.GEnergy(4.0, 'MeV'), 5.0e-5)
   spectrum.append(gammalib.GEnergy(8.0, 'MeV'), 3.0e-5)
   spectrum.save('my_file_function.txt', True)
