.. _sec_sky:

Sky maps, sky coordinates and sky regions
-----------------------------------------

Overview
~~~~~~~~

The sky module provides classes to handle sky maps, sky coordinates
and sky regions.
Sky maps are implemented by the ``GSkymap`` class which stores
intensities for a defined grid of sky coordinates.
A ``GSkymap`` object may hold several sky maps with identical
pixelisation, for example to represent the energy-dependence of the
sky intensities.
The relation between sky coordinates and pixel coordinates is
defined by the coordinate projection, represented by the abstract
``GSkyProjection`` base class. The derived class ``GHealpix``
implements the relation for the HEALPix pixelisation scheme, the 
abstract ``GWcs`` base class represents the relation for World
Coordinate Systems. Specific World Coordinate System projections
are implemented by the ``GWcsAIT``, ``GWcsAZP``, ``GWcsCAR``,
``GWcsMER``, ``GWcsSTG`` and ``GWcsTAN`` classes. Instances of
all specific World Coordinate System classes are collected in
the ``GWcsRegistry`` registry class.


Sky coordinates are implemented by the ``GSkyDir`` class that
specifies celestial coordinates in either equatorial (Right Ascension
and Declination) or galactic (longitude and latitude) coordinates.
Transformation between both systems is handled transparently by
``GSkyDir``.
Sky map pixels are implemented by the ``GSkyPixel`` class.

Sky regions are represented by the abstract ``GSkyRegion`` base class.
So far, only a simple circular sky region is implemented by the
``GSkyRegionCircle`` class. Sky regions are collected in the 
``GSkyRegions`` container class.


:ref:`fig_uml_sky` present an overview over the C++ classes of the sky
module and their relations.

.. _fig_uml_sky:

.. figure:: uml_sky.png
   :width: 100%

   Sky module
