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
#include "GTools.hpp"
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

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFits::GFits()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param fits FITS file from which the instance should be built
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


/***********************************************************************//**
 * @brief Destructor
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

/***********************************************************************//**
 * @brief  Assignment operator
 *
 * @param fits FITS file that should be assigned
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

/***********************************************************************//**
 * @brief Open FITS file
 *
 * @param filename Name of FITS file to be opened
 *
 * Opens all HDUs in the specified FITS file.
 ***************************************************************************/
void GFits::open(const std::string& filename)
{
    // Don't allow opening if another file is already open
    if (m_fitsfile != NULL)
        throw GException::fits_already_opened(G_OPEN, m_filename);

    // Open FITS file
    int status = 0;
    status     = __ffopen(&m_fitsfile, filename.c_str(), 1, &status);
    
    // If FITS file does not exist then create it now
    if (status == 104) {
        status = 0;
        status = __ffinit(&m_fitsfile, filename.c_str(), &status);
        if (status != 0)
            throw GException::fits_open_error(G_OPEN, filename, status);
    }
    else if (status != 0)
        throw GException::fits_open_error(G_OPEN, filename, status);

    // Store FITS file attributes
    m_filename = filename;

    // Determine number of HDUs
    status = __ffthdu(m_fitsfile, &m_num_hdu, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Allocate HDUs
    if (m_hdu != NULL) delete [] m_hdu;
    if (m_num_hdu > 0) m_hdu = new GFitsHDU[m_num_hdu];

    // Open all HDUs
    for (int i = 0; i < m_num_hdu; ++i)
        m_hdu[i].open(m_fitsfile, i+1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append HDU to FITS file
 *
 * @param hdu Pointer to FITS HDU that should be appended
 *
 * NOT YET FULLY IMPLEMENTED
 ***************************************************************************/
void GFits::append(const GFitsHDU* hdu)
{
    // Determine number of HDUs to add. If there are no HDUs so far and if
    // the HDU to append is not an image then we have to add a primary
    // image first. In this case we add 2 new HDUs.
    int n_add = (m_num_hdu == 0 && hdu->m_type != 0) ? 2 : 1;
    
    // Create memory to hold HDUs
    GFitsHDU* tmp = new GFitsHDU[m_num_hdu+n_add];
    if (tmp != NULL) {
    
        // Copy over existing HDUs and remove old ones
        if (m_hdu != NULL) {
            for (int i = 0; i < m_num_hdu; ++i)
                tmp[i] = m_hdu[i];
            delete [] m_hdu;
        }

        // Connect the new memory to the card pointer
        m_hdu = tmp;
        
        // Add primary image if required
        if (n_add == 2) {

            // Append empty primary image
            //m_hdu[m_num_hdu] = TBD

            // Set FITS file pointer for new HDU
            __fitsfile fptr  = *m_fitsfile;
            fptr.HDUposition = m_num_hdu;
            m_hdu[m_num_hdu].connect(&fptr);

            // Increment number of HDUs
            m_num_hdu++;
        }

        // Append new HDU to list
        m_hdu[m_num_hdu] = *hdu;
        
        // Set FITS file pointer for new HDU
        __fitsfile fptr  = *m_fitsfile;
        fptr.HDUposition = m_num_hdu;
        m_hdu[m_num_hdu].connect(&fptr);
        
        // Increment number of HDUs
        m_num_hdu++;
        
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Saves FITS file
 *
 * Saves all HDUs to the FITS file.
 ***************************************************************************/
void GFits::save(void)
{
    // Save all HDUs
    for (int i = 0; i < m_num_hdu; ++i)
        m_hdu[i].save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Saves to specified FITS file
 *
 * @param filename Name of file into which should be saved
 * @param clobber Specifies whether any existing file should be overwritten
 *
 * FUNCTION IS NOT FULLY IMPLEMENTED
 ***************************************************************************/
void GFits::saveto(const std::string& filename, int clobber)
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


/***********************************************************************//**
 * @brief Close FITS file
 *
 * Closing detaches a FITS file from the GFits object and returns a clean
 * empty object.
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


/***********************************************************************//**
 * @brief Get pointer to HDU
 *
 * @param extname Name of HDU extension which should be returned
 *
 * Returns NULL if extension has not been found.
 ***************************************************************************/
GFitsHDU* GFits::hdu(const std::string& extname)
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

/***********************************************************************//**
 * @brief Initialise class members
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


/***********************************************************************//**
 * @brief  Copy class members
 *
 * @param fits Object to be copied
 *
 * The method does not copy the FITS filename and FITS file pointer.
 * This prevents that several copies of the FITS file pointer exist in
 * different instances of GFits, which would lead to confusion since one
 * instance could close the file while for another it still would be
 * opened. The rule ONE INSTANCE - ONE FILE applies.
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


/***********************************************************************//**
 * @brief Delete class members
 *
 * This method also closes a file. In case that there are no HDUs or in case
 * that the FITS file is not correct the result file will be deleted. This
 * prevents leaving corrupted files on disk.
 ***************************************************************************/
void GFits::free_members(void)
{
    // If FITS file has been opened then close it now
    if (m_fitsfile != NULL) {
    
        // If there are no HDUs then delete the file (don't worry about error)
        if (m_num_hdu == 0) {
            int status = 0;
            __ffdelt(m_fitsfile, &status);
        }
        
        // ... otherwise close the file
        else {
            int status = 0;
            status     = __ffclos(m_fitsfile, &status);
            if (status == 252) {
                int new_status = 0;
                __ffdelt(m_fitsfile, &new_status);
                throw GException::fits_error(G_FREE_MEM, status);
            }
            else if (status != 0)
                throw GException::fits_error(G_FREE_MEM, status);
        }
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

/***********************************************************************//**
 * @brief Output operator
 *
 * @param os Output stream into which the FITS file will be dumped
 * @param fits FITS file to dump
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFits& fits)
{
    // Put header in stream
    os << "=== GFits ===" << endl;
    os << " Filename ..................: " << fits.m_filename << endl;
    os << " Number of HDUs ............: " << fits.m_num_hdu << endl;
    for (int i = 0; i < fits.m_num_hdu; ++i)
        os << fits.m_hdu[i];

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                       Other functions used by GFits                     =
 =                                                                         =
 ==========================================================================*/
