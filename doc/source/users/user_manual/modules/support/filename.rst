Filename handling
~~~~~~~~~~~~~~~~~

The :doxy:`GFilename` class handles file names in GammaLib, including 
extensions for FITS files in the cfitsio format. The following file name 
formats are supported ::

   myfile.fits
   myfile.fits[EVENTS]
   myfile.fits[EVENTS,2]
   myfile.fits[3,2]

The class decomposes the input string into the filename and the extension.
The filename can be access using the :doxy:`GFilename::filename`
method, the extension name using the :doxy:`GFilename::extname` method,
the extension number using the :doxy:`GFilename::extno` method and the 
extension version using the :doxy:`GFilename::extver` method. The latter
three methods take an argument that specifies the default value that 
should be used in case that no extension information is specified. Below
a usage example that extracts by default the ``EBOUNDS`` table from a
FITS file:

**C++**

.. code-block:: cpp
   :linenos:

   GFits     file;
   GFilename fname(filename);
   file.open(fname.filename());
   const GFitsTable& table = *file.table(fname.extname("EBOUNDS"));
