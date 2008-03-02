/***************************************************************************
 *                  GFitsImage.cpp  - FITS image class                     *
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
#include <iostream>
#include "GException.hpp"
#include "GFitsImage.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_OPEN "GFitsImage::open(fitsfile*)"
#define G_SAVE "GFitsImage::save()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                    GFitsImage constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsImage::GFitsImage() : GFitsData()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param naxis Image dimensions (0,1,2,3,4)
 * @param naxes Number of pixels in each dimension
 * @param bitpix Number of Bits per pixel (8, 16, 32, 64, -32, -64=default)
 *
 * Construct instance of GFitsImage by specifying the image dimension and
 * the number of pixels in each dimension. Note that this constructor does
 * not allocate any memory for the actual image.
 ***************************************************************************/
GFitsImage::GFitsImage(int naxis, const int* naxes, int bitpix) : GFitsData()
{
    // Initialise class members for clean destruction
    init_members();
    
    // Store number of axis and bitpix
    m_naxis  = naxis;
    m_bitpix = bitpix;
    
    // Copy number of pixels in each dimension
    if (m_naxis > 0) {
        m_naxes = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i)
            m_naxes[i] = naxes[i];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param image FITS image which should be used for construction
 ***************************************************************************/
GFitsImage::GFitsImage(const GFitsImage& image) : GFitsData(image)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(image);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsImage::~GFitsImage()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GFitsImage operators                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param image FITS image to be assigned
 ***************************************************************************/
GFitsImage& GFitsImage::operator= (const GFitsImage& image)
{
    // Execute only if object is not identical
    if (this != &image) {

        // Copy base class members
        this->GFitsData::operator=(image);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(image);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                        GFitsImage public methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Open Image
 *
 * @param fptr FITS file pointer
 *
 * Open FITS image in FITS file. If the specified image does not exist a
 * new image will be created in the FITS file.
 * IMPLEMENTATION NOT COMPLETE (IMAGE CUBE LOADING MISSING)
 ***************************************************************************/
void GFitsImage::open(__fitsfile* fptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, (fptr->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Save FITS file pointer
    m_fitsfile = *fptr;
    
    // Load image
    // TBD
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Image
 *
 * Save image pixels into FITS file. In case that the HDU does not exist it
 * is created.
 ***************************************************************************/
void GFitsImage::save(void)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
    
    // If HDU does not yet exist in file then create it now
    if (status == 107) {
        status = 0;
        status = __ffcrim(&m_fitsfile, m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE, status);
    }
    else if (status != 0)
        throw GException::fits_error(G_SAVE, status);

    // If HDU seems to be empty then create it now. This is only needed for the
    // primary HDU, since __ffmahd gives no error if the primary HDU is empty.
    // By checking the number of keywords in the HDU we detect an empty HDU ...
    int num = 0;
    status  = __ffghsp(&m_fitsfile, &num, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_SAVE, status);
    if (num == 0) {
        status = __ffcrim(&m_fitsfile, m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE, status);
    }

    // Save the image pixels (if there are some ...)
    if (m_naxis > 0) {
        long* fpixel = new long[m_naxis];
        long* lpixel = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i) {
            fpixel[i] = 1;
            lpixel[i] = m_naxes[i];
        }
        status = __ffpss(&m_fitsfile, m_type, fpixel, lpixel, m_pixels, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE, status);
        delete [] fpixel;
        delete [] lpixel;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close Image
 ***************************************************************************/
void GFitsImage::close(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        GFitsImage private methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsImage::init_members(void)
{
    // Initialise members
    m_fitsfile.HDUposition = 0;
    m_fitsfile.Fptr        = NULL;
    m_bitpix               = 8;
    m_naxis                = 0;
    m_naxes                = NULL;
    m_linked               = 0;
    m_type                 = __TBYTE;
    m_pixels               = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param image FITS image to copy
 ***************************************************************************/
void GFitsImage::copy_members(const GFitsImage& image)
{
    // Copy attributes
    m_fitsfile = image.m_fitsfile;
    m_bitpix   = image.m_bitpix;
    m_naxis    = image.m_naxis;
    m_naxes    = NULL;
    m_linked   = image.m_linked;
    m_type     = image.m_type;
    m_pixels   = NULL;

    // Copy axes
    if (image.m_naxes != NULL && m_naxis > 0) {
        m_naxes = new long[m_naxis];
        memcpy(m_naxes, image.m_naxes, m_naxis*sizeof(long));
    }
    
    // Copy pixels (not if they have been linked!!!)
    if (m_linked)
        m_pixels = image.m_pixels;
    else {
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsImage::free_members(void)
{
    // Free memory
    if (m_naxes != NULL) delete [] m_naxes;
    
    // If we have proper pixels then free them now (type dependent!)
    if (!m_linked) {
        if (m_pixels != NULL) {
            switch (m_type) {
            case __TDOUBLE:
                delete [] (double*)m_pixels;
                break;
            default:
                break;
            }
        }
    }
    
    // Mark memory as free
    m_naxes  = NULL;
    m_pixels = NULL;
    m_linked = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Connect image to FITS file
 *
 * @param fptr FITS file pointer to which the image should be connected
 ***************************************************************************/
void GFitsImage::connect(__fitsfile* fptr)
{
    // Connect Image
    m_fitsfile = *fptr;
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GFitsImage friends                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param os Output stream into which the result will be writted
 * @param image FITS image which should be put in the output stream
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitsImage& image)
{
    // Put header in stream
    os << "=== GFitsImage ===" << endl;
    os << " Number of dimensions ......: " << image.m_naxis << endl;
    for (int i = 0; i < image.m_naxis; ++i)
        os << " Number of bins in " << i << " .......: " << image.m_naxes[i] << endl;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GFitsImage                   =
 =                                                                         =
 ==========================================================================*/
