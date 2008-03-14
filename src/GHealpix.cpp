/***************************************************************************
 *             GHealpix.cpp  -  Healpix sky representation class           *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GHealpix.cpp
 * @brief GHealpix class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"
#include "GHealpix.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                     GHealpix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GHealpix::GHealpix()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param pixels GHealpix instance which should be used for construction
 ***************************************************************************/
GHealpix::GHealpix(const GHealpix& pixels)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GHealpix::~GHealpix()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GHealpix operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] pixels GHealpix instance to be assigned
 ***************************************************************************/
GHealpix& GHealpix::operator= (const GHealpix& pixels)
{
    // Execute only if object is not identical
    if (this != &pixels) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(pixels);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GHealpix public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Load Healpix data from table.
 *
 * @param hdu FITS HDU containing the Healpix data
 ***************************************************************************/
void GHealpix::load(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();
    
    // Get HDU header
    GFitsHeader* hdr = hdu->header();
    if (hdr == NULL)
        throw;

    // Check if we have a healpix representation
    if (hdr->card("PIXTYPE")->string() != "HEALPIX")
        throw;
    
    // Get Healpix resolution and determine number of pixels
    m_nside      = hdr->card("NSIDE")->integer();
    m_num_pixels = 12*m_nside*m_nside;

    // Get ordering
    try {
        std::string ordering = hdr->card("ORDERING")->string();
        if (ordering == "RING")
            m_order = 0;
        else if (ordering == "NESTED")
            m_order = 1;
        else
            throw;
    }
    catch (exception &e) {
        m_order = 1; // Default is NESTED
    }

    // Get coordinate system
    try {
        std::string coordsys = hdr->card("COORDSYS")->string();
        if (coordsys == "EQUATORIAL")
            m_coordsys = 0;
        else if (coordsys == "GALACTIC")
            m_coordsys = 1;
        else
            throw;
    }
    catch (exception &e) {
        m_coordsys = 0; // Default is EQUATORIAL
    }
    
    // Continue only of we have pixels
    if (m_num_pixels > 0) {
    
        // Get first column
        GFitsTableCol* col = hdu->column(1);

        // Check column consistency
        if (col->length() != m_num_pixels)
            throw;

        // Extract vector size of each pixel
        m_size_pixels = col->number();

        // If there are pixels then load them
        int size = m_num_pixels * m_size_pixels;
        if (size > 0) {
            m_pixels = new double[size];
            double* ptr = m_pixels;
            for (int row = 0; row < m_num_pixels; ++row) {
                for (int inx = 0; inx < m_size_pixels; ++inx, ++ptr)
                    *ptr = col->real(row,inx);
            }
        } // endif: there were pixels to load

        // Calculate longitude and latitude for each pixel
    
    } // endif: we had pixels
        
    // Return
    return;
}



/***********************************************************************//**
 * @brief Returns number of divisions of the side of each base pixel.
 ***************************************************************************/
int GHealpix::nside(void) const
{
    // Return nside
    return m_nside;
}


/***********************************************************************//**
 * @brief Returns number of pixels.
 ***************************************************************************/
int GHealpix::num_pixels(void) const
{
    // Return nside
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 ***************************************************************************/
double GHealpix::omega(void) const
{
    // Set solid angle
    double omega = fourpi / m_num_pixels;

    // Return solid angle
    return omega;
}


/*==========================================================================
 =                                                                         =
 =                         GHealpix private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GHealpix::init_members(void)
{
    // Initialise members
    m_nside       = 0;
    m_order       = 0;
    m_coordsys    = 0;
    m_num_pixels  = 0;
    m_size_pixels = 0;
    m_pixels      = NULL;
    m_lon         = NULL;
    m_lat         = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pixels GHealpix instance from which members should be copied
 ***************************************************************************/
void GHealpix::copy_members(const GHealpix& pixels)
{
    // Copy attributes
    m_nside       = pixels.m_nside;
    m_order       = pixels.m_order;
    m_coordsys    = pixels.m_coordsys;
    m_num_pixels  = pixels.m_num_pixels;
    m_size_pixels = pixels.m_size_pixels;
    
    // Copy arrays
    if (m_num_pixels > 0) {
        if (m_size_pixels > 0 && pixels.m_pixels != NULL) {
            int size = m_num_pixels*m_size_pixels;
            m_pixels = new double[size];
            memcpy(m_pixels, pixels.m_pixels, size*sizeof(double));
        }
        if (pixels.m_lon != NULL) {
            m_lon = new double[m_num_pixels];
            memcpy(m_lon, pixels.m_lon, m_num_pixels*sizeof(double));
        }
        if (pixels.m_lat != NULL) {
            m_lat = new double[m_num_pixels];
            memcpy(m_lat, pixels.m_lat, m_num_pixels*sizeof(double));
        }
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GHealpix::free_members(void)
{
    // Free memory
    if (m_pixels != NULL) delete [] m_pixels;
    if (m_lon    != NULL) delete [] m_lon;
    if (m_lat    != NULL) delete [] m_lat;

    // Mark memory as free
    m_pixels = NULL;
    m_lon    = NULL;
    m_lat    = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Healpix representation
 ***************************************************************************/
GHealpix* GHealpix::clone(void) const
{
    return new GHealpix(*this);
}


/*==========================================================================
 =                                                                         =
 =                             GHealpix friends                            =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GHealpix                    =
 =                                                                         =
 ==========================================================================*/
