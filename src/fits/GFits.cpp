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
#include <iostream>

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_OPEN     "GFits::open(std::string)"
#define G_SAVETO   "GFits::saveto(const std::string, int)"
#define G_HDU1     "GFits::hdu(const std::string&)"
#define G_HDU2     "GFits::hdu(int extno)"
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
 * @param[in] fits FITS file from which the instance should be built.
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
 * @param fits[in] FITS file that should be assigned to GFits instance.
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
 * @brief Opens or creates FITS file
 *
 * @param[in] filename Name of FITS file to be opened
 *
 * @exception GException::fits_already_opened 
 *            Class instance contains already an opened FITS file.
 *            Close file before opening a new one using GFits::close().
 * @exception GException::fits_open_error 
 *            Unable to open the specified file.
 * @exception GException::fits_error 
 *            Unable to determine number of HDUs in the FITS file.
 *
 * This method opens all HDUs that are found in the specified FITS file.
 * If the file does not exist then a new FITS file is created.
 * For each HDU, a GFitsHDU object is associated to the GFits object.
 * The HDUs can then be accessed using the
 * GFits::hdu(const std::string&)
 * or
 * GFits::hdu(int extno)
 * methods.
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
 * @param[in] hdu FITS HDU that should be appended
 *
 * Append HDU to the next free position in a FITS file. In case that no HDU
 * exists so far in the FITS file and if the HDU to append is NOT an image,
 * an empty primary image will be inserted as first HDU in the FITS file.
 * This guarantees the compatibility with the FITS standard.
 ***************************************************************************/
void GFits::append_hdu(const GFitsHDU& hdu)
{
    // Determine number of HDUs to add. If there are no HDUs so far and if
    // the HDU to append is not an image then we have to add a primary
    // image first. In this case we add 2 new HDUs.
    int n_add = (m_num_hdu == 0 && hdu.m_type != 0) ? 2 : 1;

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
            m_hdu[m_num_hdu].primary();

            // Set FITS file pointer for new HDU
            __fitsfile fptr  = *m_fitsfile;
            fptr.HDUposition = m_num_hdu;
            m_hdu[m_num_hdu].connect(&fptr);

            // Increment number of HDUs
            m_num_hdu++;
        }

        // Append new HDU to list
        m_hdu[m_num_hdu] = hdu;

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
 * Saves all HDUs to the FITS file by looping over the GFitsHDU::save()
 * method.
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
 * @param[in] filename Name of file into which should be saved.
 * @param[in] clobber Specifies whether any existing file should be 
 *            overwritten.
 *
 * @exception GException::fits_file_exist
 *            File specified by 'filename' exists already.
 *            To overwrite an existing file set 'clobber=1'.
 * @exception GException::fits_error
 *            Unable to delete exiting FITS file.
 ***************************************************************************/
void GFits::saveto(const std::string& filename, int clobber)
{
    // Create a copy of existing object. Recall that the copy does not
    // carry over the FITS filename and pointer. Those will be automatically
    // reset to the initial values in the new object.
    GFits new_fits = *this;

    // Check if specified FITS file exists. If yes, saving will only be
    // allowed if clobber is true. If this is the case the specified file
    // will be deleted.
    int status = 0;
    status     = __ffopen(&(new_fits.m_fitsfile), filename.c_str(), 1, &status);
    if (status == 0) {

        // If overwriting was not allowed the throw an exception
        if (!clobber)
            throw GException::fits_file_exist(G_SAVETO, filename, status);

        // Delete existing file
        __fitsfile* tmp;
        status = __ffopen(&tmp, filename.c_str(), 1, &status);
        status = __ffdelt(tmp, &status);
        if (status == 0)
            throw GException::fits_error(G_SAVETO, status);

    }
    else
        status = 0;

    // Create a new FITS file now
    new_fits.open(filename);

    // Copy all HDUs
    for (int i = 0; i < m_num_hdu; ++i)
        new_fits.m_hdu[i] = m_hdu[i];

    // Save new FITS file
    new_fits.save();

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
 * @param[in] extname Name of HDU extension which should be returned
 *
 * @exception GException::fits_hdu_not_found
 *            HDU with specified extension name has not been found in FITS
 *            file.
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

    // Throw an error if HDU has not been found
    if (ptr == NULL)
        throw GException::fits_hdu_not_found(G_HDU1, extname);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Get pointer to HDU
 *
 * @param extno[in] Extension number (starting from 1)
 *
 * @exception GException::fits_hdu_not_found
 *            HDU with specified extension number has not been found in FITS
 *            file.
 ***************************************************************************/
GFitsHDU* GFits::hdu(int extno) 
{
    // Verify that extension number is in valid range
    if (extno < 1 || extno > m_num_hdu)
        throw GException::out_of_range("GFits::hdu(int)", extno, 1, m_num_hdu);

    // Get HDU pointer
    GFitsHDU* ptr = &(m_hdu[extno-1]);

    // Throw an error if HDU has not been found
    if (ptr == NULL) {
        ostringstream s_extname;
        s_extname << extno;
        std::string extname = "extno=" + s_extname.str();
        throw GException::fits_hdu_not_found(G_HDU2, extname);
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
