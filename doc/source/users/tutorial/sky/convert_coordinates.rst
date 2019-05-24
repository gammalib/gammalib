How to convert a sky coordinate from celestial to Galactic?
===========================================================

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
