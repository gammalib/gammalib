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
#define G_OPEN     "GFits::open(std::string)"
#define G_SAVETO   "GFits::saveto(const std::string, int)"
#define G_FREE_MEM "GFits::free_members()"

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
void GFits::open(const std::string filename)
{
    // Don't allow opening if another file is already open
    if (m_fitsfile != NULL)
        throw GException::fits_already_opened(G_OPEN, m_filename);

    // Open FITS file
    int status = 0;
    status     = __ffopen(&m_fitsfile, filename.c_str(), 1, &status);
    if (status != 0)
        throw GException::fits_open_error(G_OPEN, filename, status);

    // Store FITS file attributes
    m_filename = filename;

    // Determine number of HDUs
    status = __ffthdu(m_fitsfile, &m_num_hdu, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

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
 *                               Save FITS file                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::save(void)
{
    // Save all HDUs
    for (int i = 0; i < m_num_hdu; ++i)
        m_hdu[i].save();

    // Return
    return;
}


/***************************************************************************
 *                       Save to specified FITS file                       *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFits::saveto(const std::string filename, int clobber)
{
    // Create a copy of existing object. Recall that the copy does not
    // carry over the FITS filename and pointer. Those will be automatically
    // reset to the initial values in the new object.
    GFits new_fits = *this;

    // Check if specified FITS file exists. If yes, saving will only be
    // allowed if clobber is true ...
    int status = 0;
    status     = __ffopen(&(new_fits.m_fitsfile), filename.c_str(), 1, &status);
    if (status == 0) {
        if (!clobber) {
            throw GException::fits_open_error(G_SAVETO, filename, status);
            //throw file_exists ...
        }
        //delete existing file
        //...
    }

    // ... otherwise create a new FITS file now
    else {
        // ...
    }

    // Save all HDUs
    for (int i = 0; i < m_num_hdu; ++i)
        new_fits.m_hdu[i].save();

    // Return
    return;
}


/***************************************************************************
 *                              Close FITS file                            *
 * ----------------------------------------------------------------------- *
 * Closing detaches a FITS file from the GFits object and returns a clean  *
 * empty object.                                                           *
 ***************************************************************************/
void GFits::close(void)
{
    // Close file and free all members 
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***************************************************************************
 *                            Get pointer to HDU                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsHDU* GFits::hdu(std::string extname)
{
    // Initialise result to NULL pointer
    GFitsHDU* ptr = NULL;

    // Search for specified extension
    for (int i = 0; i < m_num_hdu; ++i) {
        if (m_hdu[i].extname() == extname) {
            ptr = &(m_hdu[i]);
            break;
        }
    }

    // Return pointer
    return ptr;
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
 * The function does not copy the FITS filename and FITS file pointer.     *
 * This prevents that several copies of the FITS file pointer exist in     *
 * different instances of GFits, which would lead to confusion since one   *
 * instance could close the file while for another it still would be       *
 * opened. The rule ONE INSTANCE - ONE FILE applies.                       *
 ***************************************************************************/
void GFits::copy_members(const GFits& fits)
{
    // Reset FITS file attributes
    m_filename.clear();
    m_fitsfile = NULL;

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
    // If FITS file has been opened then close it now
    if (m_fitsfile != NULL) {
        int status = 0;
        status     = __ffclos(m_fitsfile, &status);
        if (status != 0)
            throw GException::fits_error(G_FREE_MEM, status);
    }

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
