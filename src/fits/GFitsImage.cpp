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

/* __ Method name definitions ____________________________________________ */
#define G_NAXES      "GFitsImage::naxes(int)"
#define G_OPEN_IMAGE "GFitsImage::open(fitsfile*)"
#define G_LOAD_IMAGE "GFitsImage::load_image(int, const void*, const void*, int*)"
#define G_SAVE_IMAGE "GFitsImage::save_image(int, const void*)"

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
 *
 * Construct instance of GFitsImage by specifying the image dimension and
 * the number of pixels in each dimension. Note that this constructor does
 * not allocate any memory for the actual image.
 ***************************************************************************/
GFitsImage::GFitsImage(int naxis, const int* naxes) : GFitsData()
{
    // Initialise class members for clean destruction
    init_members();
    
    // Store number of axis and bitpix
    m_naxis = naxis;
    
    // Copy number of pixels in each dimension and calculate the total
    // number of pixels
    if (m_naxis > 0) {
        m_naxes      = new long[m_naxis];
        m_num_pixels = 1;
        for (int i = 0; i < m_naxis; ++i) {
            m_naxes[i]    = naxes[i];
            m_num_pixels *= naxes[i];
        }
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
 * @param[in] image FITS image to be assigned
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
 * @brief Return number of Bits per pixel (negative=floating point)
 ***************************************************************************/
int GFitsImage::bitpix(void) const
{
    // Return number of Bits per pixel
    return m_bitpix;
}


/***********************************************************************//**
 * @brief Return image dimension
 ***************************************************************************/
int GFitsImage::naxis(void) const
{
    // Return dimension
    return m_naxis;
}


/***********************************************************************//**
 * @brief Return axis dimension
 *
 * @param[in] axis Axis for which the dimension should be returned 
 *            (starting from 0)
 ***************************************************************************/
int GFitsImage::naxes(int axis) const
{
    // Check if axis is within the range
    if (axis < 0 || axis >= m_naxis)
        throw GException::out_of_range(G_NAXES, axis, 0, m_naxis-1);

    // Get axis dimension
    int dim = m_naxes[axis];

    // Return axis dimension
    return dim;
}


/***********************************************************************//**
 * @brief Return number of pixels
 ***************************************************************************/
int GFitsImage::num_pixels(void) const
{
    // Return number of pixels
    return m_num_pixels;
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
    m_num_pixels           = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] image FITS image to copy
 ***************************************************************************/
void GFitsImage::copy_members(const GFitsImage& image)
{
    // Copy attributes
    m_fitsfile   = image.m_fitsfile;
    m_bitpix     = image.m_bitpix;
    m_naxis      = image.m_naxis;
    m_naxes      = NULL;
    m_num_pixels = image.m_num_pixels;

    // Copy axes
    if (image.m_naxes != NULL && m_naxis > 0) {
        m_naxes = new long[m_naxis];
        memcpy(m_naxes, image.m_naxes, m_naxis*sizeof(long));
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

    // Mark memory as free
    m_naxes = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Open Image
 *
 * @param[in] fptr FITS file pointer
 *
 * Open FITS image in FITS file. Opening means connecting the FITS file
 * pointer to the image and reading the image and axes dimensions.
 ***************************************************************************/
void GFitsImage::open_image(__fitsfile* fptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, (fptr->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN_IMAGE, status);

    // Save FITS file pointer
    m_fitsfile = *fptr;

    // Get the image dimensions
    status = __ffgidm(&m_fitsfile, &m_naxis, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN_IMAGE, status);

    // Reset number of image pixels
    m_num_pixels = 0;

    // Get the axes dimensions
    if (m_naxis > 0) {

        // Allocate memory for axes dimensions
        if (m_naxes != NULL) delete [] m_naxes;
        m_naxes      = new long[m_naxis];

        // Get the axes dimensions
        status = __ffgisz(&m_fitsfile, m_naxis, m_naxes, &status);
        if (status != 0)
            throw GException::fits_error(G_OPEN_IMAGE, status);

        // Calculate number of image pixels
        m_num_pixels = 1;
        for (int i = 0; i < m_naxis; ++i)
            m_num_pixels *= m_naxes[i];

    } // endif: there is an image

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load FITS image
 *
 * @param[in] datatype Datatype of pixels to be saved
 * @param[in] pixels Pixel array to be saved
 *
 * Save image pixels into FITS file. In case that the HDU does not exist it
 * is created.
 ***************************************************************************/
void GFitsImage::load_image(int datatype, const void* pixels, const void* nulval,
                            int* anynul)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_LOAD_IMAGE, status);

    // Load the image pixels (if there are some ...)
    if (m_naxis > 0) {
        long* fpixel = new long[m_naxis];
        long* lpixel = new long[m_naxis];
        long* inc    = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i) {
            fpixel[i] = 1;
            lpixel[i] = m_naxes[i];
            inc[i]    = 1;
        }
        status = __ffgsv(&m_fitsfile, datatype, fpixel, lpixel, inc, (void*)nulval,
                         (void*)pixels, anynul, &status);
        if (status != 0)
            throw GException::fits_error(G_LOAD_IMAGE, status);
        delete [] fpixel;
        delete [] lpixel;
        delete [] inc;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save FITS image
 *
 * @param[in] datatype Datatype of pixels to be saved
 * @param[in] pixels Pixel array to be saved
 *
 * Save image pixels into FITS file. In case that the HDU does not exist it
 * is created. In case that the pixel array is empty no data are saved; all
 * image pixels will be empty in this case.
 ***************************************************************************/
void GFitsImage::save_image(int datatype, const void* pixels)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);

    // If HDU does not yet exist in file then create it now
    if (status == 107) {
        status = 0;
        status = __ffcrim(&m_fitsfile, m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE_IMAGE, status);
    }
    else if (status != 0)
        throw GException::fits_error(G_SAVE_IMAGE, status);

    // If HDU seems to be empty then create it now. This is only needed for the
    // primary HDU, since __ffmahd gives no error if the primary HDU is empty.
    // By checking the number of keywords in the HDU we detect an empty HDU ...
    int num = 0;
    status  = __ffghsp(&m_fitsfile, &num, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_SAVE_IMAGE, status);
    if (num == 0) {
        status = __ffcrim(&m_fitsfile, m_bitpix, m_naxis, m_naxes, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE_IMAGE, status);
    }

    // Save the image pixels (if there are some ...)
    if (m_naxis > 0 && pixels != NULL) {
        long* fpixel = new long[m_naxis];
        long* lpixel = new long[m_naxis];
        for (int i = 0; i < m_naxis; ++i) {
            fpixel[i] = 1;
            lpixel[i] = m_naxes[i];
        }
        status = __ffpss(&m_fitsfile, datatype, fpixel, lpixel, (void*)pixels, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE_IMAGE, status);
        delete [] fpixel;
        delete [] lpixel;
    }

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
std::ostream& operator<< (std::ostream& os, const GFitsImage& image)
{
    // Put header in stream
    os << "=== GFitsImage ===" << std::endl;
    os << " Number of dimensions ......: " << image.m_naxis << std::endl;
    os << " Number of image pixels ....: " << image.m_num_pixels << std::endl;
    for (int i = 0; i < image.m_naxis; ++i)
        os << " Number of bins in " << i << " .......: " << image.m_naxes[i] << std::endl;

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                    Other functions used by GFitsImage                   =
 =                                                                         =
 ==========================================================================*/
