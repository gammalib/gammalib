/***************************************************************************
 *             GSkymap.cpp  -  Class that implements a sky map             *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2010 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
//#include <cmath>
#include "GException.hpp"
//#include "GTools.hpp"
#include "GSkymap.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototype __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                      GSkymap constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GSkymap::GSkymap(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix constructor
 *
 * @param[in] wcs World Coordinate System (HPX).
 * @param[in] coordsys Coordinate System (CEL or GAL).
 * @param[in] nside Nside parameter.
 * @param[in] ordering Pixel ordering (RING or NEST).
 * @param[in] nmaps Number of maps in set (default=1).
 ***************************************************************************/
GSkymap::GSkymap(const std::string& wcs, const std::string& coordsys, 
                 const int& nside, const std::string& ordering,
                 const int nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    // TODO

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix constructor
 *
 * @param[in] wcs World Coordinate System (HPX).
 * @param[in] coordsys Coordinate System (CEL or GAL).
 * @param[in] dir Centre of skymap.
 * @param[in] nlon Number of pixels in longitude.
 * @param[in] nlat Number of pixels in latitude.
 * @param[in] dlon Size of pixels in longitude.
 * @param[in] dlat Size of pixels in latitude.
 * @param[in] nmaps Number of maps in set (default=1).
 ***************************************************************************/
GSkymap::GSkymap(const std::string& wcs, const std::string& coordsys, 
                 GSkyDir& dir, const int& nlon, const int& nlat,
                 const double& dlon, const double& dlat, const int nmaps)
{
    // Initialise class members for clean destruction
    init_members();
    
    // TODO

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] map Sky map from which class should be instantiated.
 ***************************************************************************/
GSkymap::GSkymap(const GSkymap& map)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(map);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkymap::~GSkymap(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GSkymap operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 ***************************************************************************/
GSkymap& GSkymap::operator= (const GSkymap& map)
{
    // Execute only if object is not identical
    if (this != &map) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(map);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GSkymap public methods                         =
 =                                                                         =
 ==========================================================================*/



/*==========================================================================
 =                                                                         =
 =                          GSkymap private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkymap::init_members(void)
{
    // Initialise members
    m_coordsys   = 0;
    m_num_pixels = 0;
    m_wcs        = NULL;
    m_pixels     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Sky map from which members should be copied
 ***************************************************************************/
void GSkymap::copy_members(const GSkymap& map)
{
    // Copy attributes
    m_coordsys   = map.m_coordsys;
    m_num_pixels = map.m_num_pixels;
    
    // Copy WCS
    m_wcs = map.m_wcs->clone();
    
    // Copy pixels
    if (m_num_pixels > 0 && map.m_pixels != NULL) {
        m_pixels = new double[m_num_pixels];
        for (int i = 0; i <  m_num_pixels; ++i)
            m_pixels[i] = map.m_pixels[i];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkymap::free_members(void)
{
    // Free memory
    if (m_pixels != NULL) delete [] m_pixels;

    // Signal free pointers
    m_pixels     = NULL;
    m_num_pixels = 0;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GSkymap friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] map Sky map to put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSkymap& map)
{
    // Put header in stream
    os << "=== GSkymap ===" << std::endl;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GSkymap                    =
 =                                                                         =
 ==========================================================================*/
