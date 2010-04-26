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
#include "GFitsTableDblCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                               "GSkymap::read(const GFitsHDU*)"

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

/***********************************************************************//**
 * @brief Read skymap from FITS table.
 *
 * @param[in] hdu FITS HDU containing the Healpix data
 ***************************************************************************/
void GSkymap::read(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();
    
    // TODO: Detect here the skymap type (table = Healpix)
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read skymap from FITS table.
 *
 * @param[in] fits FITS file into which the skymap will be written
 ***************************************************************************/
void GSkymap::write(const GFits* file)
{
    // TODO: Skymap type writing


    // Append HDU to FITS file
//    fits->append_hdu(hdu);

    // Return
    return;
}


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
    m_num_pixels = 0;
    m_num_maps   = 0;
    m_wcs        = NULL;
    m_pixels     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate class members
 ***************************************************************************/
void GSkymap::alloc_pixels(void)
{
    // Compute data size
    int size = m_num_pixels * m_num_maps;

    // Continue only if there are pixels    
    if (size > 0) {
    
        // Allocate pixels and initialize them to 0
        m_pixels = new double[size];
        for (int i = 0; i < size; ++i)
            m_pixels[i] = 0.0;

    } // endif: there were pixels
    
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
    m_num_pixels = map.m_num_pixels;
    m_num_maps   = map.m_num_maps;
    
    // Copy WCS
    m_wcs = map.m_wcs->clone();

    // Compute data size
    int size = m_num_pixels * m_num_maps;
    
    // Copy pixels
    if (size > 0) {
        alloc_pixels();
        for (int i = 0; i <  size; ++i)
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
    if (m_wcs    != NULL) delete m_wcs;
    if (m_pixels != NULL) delete [] m_pixels;

    // Signal free pointers
    m_wcs        = NULL;
    m_pixels     = NULL;
    m_num_pixels = 0;
    m_num_maps   = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Healpix data from FITS table.
 *
 * @param[in] hdu FITS HDU containing the Healpix data
 ***************************************************************************/
void GSkymap::read_healpix(const GFitsHDU* hdu)
{
    // Allocate Healpix WCS
    m_wcs = new GWcsHPX;
    
    // Read WCS information from FITS header
    m_wcs->read(hdu);

    // Get first column
    GFitsTableCol* col = hdu->column(0);

    // Extract pixel information
    m_num_pixels = col->length();
    m_num_maps   = col->number();
    
    // If we have pixels then read them now
    int size = m_num_pixels * m_num_maps;
    if (size > 0) {

        // Allocate pixels
        alloc_pixels();

        // Read pixels
        double* ptr = m_pixels;
        for (int row = 0; row < m_num_pixels; ++row) {
            for (int inx = 0; inx < m_num_maps; ++inx, ++ptr)
                *ptr = col->real(row,inx);
        }
    
    } // endif: we had pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create FITS HDU containing Healpix data
 ***************************************************************************/
GFitsHDU* GSkymap::create_hdu_healpix(void)
{
    // Initialise result to NULL pointer
    GFitsHDU* hdu = NULL;
    
    // Compute size of Healpix data
    int size = m_num_pixels * m_num_maps;

    // Continue only if we have a WCS and some data
    if (m_wcs != NULL && size > 0) {
    
        // Create binary table with one column
        GFitsBinTable    table  = GFitsBinTable(m_num_pixels);
        GFitsTableDblCol column = GFitsTableDblCol("DATA", m_num_pixels, 
                                                    m_num_maps);  

        // Fill data into column
        double* ptr = m_pixels;
        for (int row = 0; row < m_num_pixels; ++row) {
            for (int inx = 0; inx < m_num_maps; ++inx, ++ptr)
                column(row,inx) = *ptr;
        }

        // Append column to table
        table.append_column(column);

        // Create HDU with table
        hdu = new GFitsHDU(table);

        // Set extension name
        hdu->extname("HEALPIX");

        // Write WCS information into FITS header
        m_wcs->write(hdu);
    
    }

    // Return HDU
    return hdu;
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
