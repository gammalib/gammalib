/***************************************************************************
 *                    GFits.cpp  - FITS file access class                  *
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

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GFits.hpp"
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                       GFits constructors/destructors                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFits::GFits()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***************************************************************************
 *                              Copy constructor                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFits::GFits(const GFits& fits)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(fits);

    // Return
    return;
}


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFits::~GFits()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GFits operators                             =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFits& GFits::operator= (const GFits& fits)
{
    // Execute only if object is not identical
    if (this != &fits) {
  
        // Free members
        free_members();
  
        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(fits);
	
    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GFits public methods                           =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                              Open FITS file                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::open(std::string filename)
{
    // Open FITS file
    int status = 0;
    status     = __ffopen(&m_fitsfile, filename.c_str(), 1, &status);
    if (status != 0) {
        throw GException::fits_open_error("GFits::open(std::string)", filename, status);
    }

    // Store FITS file attributes
    m_filename = filename;
    
    // Determine number of HDUs
    status = __ffthdu(m_fitsfile, &m_num_hdu, &status);
    if (status != 0) {
        throw GException::fits_error("GFits::open(std::string)", status);
    }
    
    // Allocate HDUs
    if (m_hdu != NULL) delete [] m_hdu;
    m_hdu = new GFitsHDU[m_num_hdu];

    // Open all HDUs
    for (int i = 0; i < m_num_hdu; ++i)
        m_hdu[i].open(m_fitsfile, i+1);
    
    // Return
    return;
}


/***************************************************************************
 *                              Close FITS file                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::close(void)
{
    // Close FITS file
    int status = 0;
    status     = __ffclos(m_fitsfile, &status);
    
    // Handle error
    if (status != 0) {
        throw GException::fits_error("GFits::close(std::string)", status);
    }
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFits private methods                          =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::init_members(void)
{
    // Initialise GFits members
    m_filename.clear();
    m_fitsfile = NULL;
    m_num_hdu  = 0;
    m_hdu      = NULL;
  
    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::copy_members(const GFits& fits)
{
    // Copy GFits attributes
    m_filename = fits.m_filename;
    m_fitsfile = fits.m_fitsfile;
    
    // Copy HDUs
    if (fits.m_hdu != NULL && fits.m_num_hdu > 0) {
        m_num_hdu = fits.m_num_hdu;
        m_hdu     = new GFitsHDU[fits.m_num_hdu];
        for (int i = 0; i < fits.m_num_hdu; ++i)
            m_hdu[i] = fits.m_hdu[i];
    }
    
    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::free_members(void)
{
    // Free memory
    if (m_hdu != NULL) delete [] m_hdu;

    // Properly mark as free
    m_hdu = NULL;
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               GFits friends                             =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GFits                     =
 =                                                                         =
 ==========================================================================*/
