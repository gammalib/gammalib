Sky coordinates, sky maps and sky regions
-----------------------------------------

How to convert a sky coordinate from celestial to Galactic?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To convert from celestial to Galactic coordinates you can set a sky 
direction in celestial coordinates and read it back in Galactic
coordinates.

Python:

.. code-block:: python

    >>> import gammalib
    >>> dir=gammalib.GSkyDir()
    >>> dir.radec_deg(83.63,22.01)
    >>> l=dir.l_deg()
    >>> b=dir.b_deg()
    >>> print(l,b)
    (184.55973405309402, -5.7891829467816827)

C++:

.. code-block:: cpp

    #include "GammaLib.hpp"
    GSkyDir dir;
    dir.radec_deg(83.63,22.01);
    double l = dir.l_deg();
    double b = dir.b_deg();
    std::cout << l << ", " << b << std::endl;


How to convert sky map projections?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following code illustrates how to convert a map in HealPix projection
into a map in cartesian projection. Map conversion is performed using the
+= operator that adds the bilinearly interpolated intensity values from one
map to another. The example code applies to any kind of map projections.

Python:

.. code-block:: python

    >>> import gammalib
    >>> healpix = gammalib.GSkyMap("healpix.fits")
    >>> map = gammalib.GSkyMap("CAR","GAL",0.0,0.0,0.5,0.5,100,100)
    >>> map += healpix
    >>> map.save("carmap.fits")

C++:

.. code-block:: cpp

    #include "GammaLib.hpp"
    GSkyMap healpix("healpix.fits");
    GSkyMap map("CAR","GAL",0.0,0.0,0.5,0.5,100,100);
    map += healpix;
    map.save("carmap.fits");
    