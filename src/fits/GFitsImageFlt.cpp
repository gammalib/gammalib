/***************************************************************************
 *         GFitsImageFlt.cpp  - FITS single precision image class          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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
#include "GFitsImageFlt.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR2 "GFitsImageFlt::operator() (int, int)"
#define G_OPEN      "GFitsImageFlt::open(fitsfile*)"
#define G_SAVE      "GFitsImageFlt::save()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                  GFitsImageFlt constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsImageFlt::GFitsImageFlt() : GFitsImage()
{
    // Initialise class members for clean destruction
    init_members();

    // Set bitpix
    m_bitpix = -32;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param naxis Image dimension (0,1,2,3,4)
 * @param naxes Number of pixels in each dimension
 *
 * Construct instance of GFitsImageFlt by specifying the image dimension and
 * the number of pixels in each dimension. Note that this constructor does
 * not allocate any memory for the actual image.
 ***************************************************************************/
GFitsImageFlt::GFitsImageFlt(int naxis, const int* naxes) : GFitsImage(naxis, naxes)
{
    // Initialise class members for clean destruction
    init_members();

    // Set bitpix
    m_bitpix = -32;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param naxis Image dimension (0,1,2,3,4)
 * @param naxes Number of pixels in each dimension
 * @param pixels Pointer to image pixel array
 *
 * Construct instance of GFitsImageFlt by specifying the image dimension, the
 * number of pixels in each dimension and an image pixel array. If the
 * 'pixels' pointer is valid (i.e. not NULL) the pixel data is indeed copied
 * by this constructor. In case of memory shortage a combination of
 *   GFitsImageFlt::GFitsImageFlt(int naxis, const int* naxes)
 * and
 *   GFitsImageFlt::link(double* pixels)
 * may be prefered to avoid duplication of the pixel data.
 ***************************************************************************/
GFitsImageFlt::GFitsImageFlt(int naxis, const int* naxes, const float* pixels) 
                                                   : GFitsImage(naxis, naxes)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy pixels into object
    if (m_num_pixels > 0 && pixels != NULL) {
        m_pixels = new float[m_num_pixels];
        for (int i = 0; i < m_num_pixels; ++i)
            m_pixels[i] = pixels[i];
    }

    // Set bitpix
    m_bitpix = -32;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param image FITS image which should be used for construction
 ***************************************************************************/
GFitsImageFlt::GFitsImageFlt(const GFitsImageFlt& image) : GFitsImage(image)
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
GFitsImageFlt::~GFitsImageFlt()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsImageFlt operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param image FITS image to be assigned
 ***************************************************************************/
GFitsImageFlt& GFitsImageFlt::operator= (const GFitsImageFlt& image)
{
    // Execute only if object is not identical
    if (this != &image) {

        // Copy base class members
        this->GFitsImage::operator=(image);

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


/***********************************************************************//**
 * @brief 1D image pixel access operator
 *
 * @param[in] ix Pixel index
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 1D image
 *
 * Provides access to a pixel of a 1D image. No range checking is performed.
 * Use the at(ix) method if range checking is required.
 ***************************************************************************/
float& GFitsImageFlt::operator() (const int& ix)
{
    // Operator is only valid for 1D images
    if (m_naxis != 1)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 1);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_pixels();

    // Return image pixel
    return m_pixels[ix];
}


/***********************************************************************//**
 * @brief 1D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 1D image
 *
 * Provides access to a pixel of a 1D image. No range checking is performed.
 * Use the at(ix) method if range checking is required.
 ***************************************************************************/
const float& GFitsImageFlt::operator() (const int& ix) const
{
    // Operator is only valid for 1D images
    if (m_naxis != 1)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 1);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageFlt*)this)->fetch_pixels();

    // Return image pixel
    return m_pixels[ix];
}


/***********************************************************************//**
 * @brief 2D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension
 * @param[in] iy Pixel index in second dimension
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 2D image
 *
 * Provides access to a pixel of a 2D image. No range checking is performed.
 * Use the at(ix,iy) method if range checking is required.
 ***************************************************************************/
float& GFitsImageFlt::operator() (const int& ix, const int& iy)
{
    // Operator is only valid for 2D images
    if (m_naxis != 2)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 2);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_pixels();

    // Calculate pixel offset
    int offset = ix + iy * m_naxes[0];

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 2D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension
 * @param[in] iy Pixel index in second dimension
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 2D image
 *
 * Provides access to a pixel of a 2D image. No range checking is performed.
 * Use the at(ix,iy) method if range checking is required.
 ***************************************************************************/
const float& GFitsImageFlt::operator() (const int& ix, const int& iy) const
{
    // Operator is only valid for 2D images
    if (m_naxis != 2)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 2);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageFlt*)this)->fetch_pixels();

    // Calculate pixel offset
    int offset = ix + iy * m_naxes[0];

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 3D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension
 * @param[in] iy Pixel index in second dimension
 * @param[in] iz Pixel index in third dimension
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 3D image
 *
 * Provides access to a pixel of a 3D image. No range checking is performed.
 * Use the at(ix,iy,iz) method if range checking is required.
 ***************************************************************************/
float& GFitsImageFlt::operator() (const int& ix, const int& iy, const int& iz)
{
    // Operator is only valid for 3D images
    if (m_naxis != 3)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 3);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_pixels();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + iz * m_naxes[1]);

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 3D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension
 * @param[in] iy Pixel index in second dimension
 * @param[in] iz Pixel index in third dimension
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 3D image
 *
 * Provides access to a pixel of a 3D image. No range checking is performed.
 * Use the at(ix,iy,iz) method if range checking is required.
 ***************************************************************************/
const float& GFitsImageFlt::operator() (const int& ix, const int& iy, const int& iz) const
{
    // Operator is only valid for 3D images
    if (m_naxis != 3)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 3);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageFlt*)this)->fetch_pixels();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + iz * m_naxes[1]);

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 4D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension
 * @param[in] iy Pixel index in second dimension
 * @param[in] iz Pixel index in third dimension
 * @param[in] it Pixel index in forth dimension
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 4D image
 *
 * Provides access to a pixel of a 4D image. No range checking is performed.
 * Use the at(ix,iy,iz,it) method if range checking is required.
 ***************************************************************************/
float& GFitsImageFlt::operator() (const int& ix, const int& iy, const int& iz, const int& it)
{
    // Operator is only valid for 4D images
    if (m_naxis != 4)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 4);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_pixels();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2]));

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 4D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension
 * @param[in] iy Pixel index in second dimension
 * @param[in] iz Pixel index in third dimension
 * @param[in] it Pixel index in forth dimension
 *
 * @exception GException::fits_wrong_image_operator
 *            Image is not a 4D image
 *
 * Provides access to a pixel of a 4D image. No range checking is performed.
 * Use the at(ix,iy,iz,it) method if range checking is required.
 ***************************************************************************/
const float& GFitsImageFlt::operator() (const int& ix, const int& iy, const int& iz, const int& it) const
{
    // Operator is only valid for 4D images
    if (m_naxis != 4)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 4);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageFlt*)this)->fetch_pixels();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2]));

    // Return image pixel
    return m_pixels[offset];
}


/*==========================================================================
 =                                                                         =
 =                       GFitsImageFlt public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Link external pixel array to image
 *
 * @param pixels Pointer to external pixel array that should be linked
 *
 * This method allows to use an externally allocated pixel array within
 * the class. No control over the existence of this pixel array is performed
 * hence the user has to assure its existence over the entire lifetime of
 * the class instance. Using of an externally allocated array may be
 * preferred in case of large image arrays that are created outside the
 * class. In this case linking prevents a duplication of the image array.
 * In all other cases, linking should be better avoided. 
 ***************************************************************************/
void GFitsImageFlt::link(float* pixels)
{
    // Free any existing pixels
    if (!m_linked) {
        if (m_pixels != NULL) delete [] m_pixels;
    }

    // Link pixel array to image
    m_pixels = pixels;

    // Flag that pixels have been linked
    m_linked = 1;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set NULL value
 *
 * @param value Value to which NULL pixels should be set
 *
 * NOTE: FOR THE MOMENT THE IMAGE IS NOT RE-LOADED IF THE NULL VALUE IS
 * CHANGED !!!!
 ***************************************************************************/
void GFitsImageFlt::set_nullval(const float* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new float;
        *m_nulval = *value;
    }

    // Re-load image 
    // NOTE: THIS WILL LEAD TO A LOSS OF MODIFICATIONS; ISSUE SAVE BEFORE !!!
    //if (m_data != NULL) {
    //    //save();
    //    load();
    //}

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer to image pixel
 ***************************************************************************/
float* GFitsImageFlt::pixels(void)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_pixels();

    // Return
    return m_pixels;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsImageFlt private methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsImageFlt::init_members(void)
{
    // Initialise members
    m_linked = 0;
    m_pixels = NULL;
    m_anynul = 0;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param image FITS image to copy
 *
 * Copy class members. If the pixel array is linked the pointer to the array
 * is copied. Otherwise a copy of the image pixels will be created.
 ***************************************************************************/
void GFitsImageFlt::copy_members(const GFitsImageFlt& image)
{
    // Copy attributes
    m_linked = image.m_linked;
    m_anynul = image.m_anynul;

    // Copy pixels
    if (m_linked)
        m_pixels = image.m_pixels;
    else {
        if (m_num_pixels > 0 && image.m_pixels != NULL) {
            m_pixels = new float[m_num_pixels];
            for (int i = 0; i < m_num_pixels; ++i)
                m_pixels[i] = image.m_pixels[i];
        }
    }

    // Copy NULL value
    if (image.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new float;
        *m_nulval = *image.m_nulval;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsImageFlt::free_members(void)
{
    // If we have proper pixels then free them now
    if (!m_linked) {
        if (m_pixels != NULL) delete [] m_pixels;
    }

    // Free NULL value
    if (m_nulval != NULL) delete m_nulval;

    // Mark memory as free
    m_linked = 0;
    m_pixels = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch image pixels
 *
 * Fetch the image pixels. This function is in general called if pixel
 * values should be read or written yet no pixel array is allocated. In case
 * that pixels existed already before they will be deleted before fetching
 * new ones. 
 * There are two possibilities to fetch the pixels:
 * (1) In case that a FITS file is attached to the image the pixel array
 * will be loaded from the file.
 * (2) In case that no FITS file is attached a new pixel array will be
 * allocated that is initalised to 0.
 ***************************************************************************/
void GFitsImageFlt::fetch_pixels(void)
{
    // Fetch only if there are pixels in image
    if (m_num_pixels > 0) {

        // If we have pixels allocated then free them now
        if (!m_linked) {
            if (m_pixels != NULL) delete [] m_pixels;
        }

        // Allocate fresh image pixels
        m_pixels = new float[m_num_pixels];

        // Case 1: A FITS is attached. Load pixels now.
        if (m_fitsfile.Fptr != NULL)
            load_image(__TFLOAT, m_pixels, m_nulval, &m_anynul);

        // Case 2: No FITS file is attached. Allocate pixels and set them all to 0.0
        else {
            for (int i = 0; i < m_num_pixels; ++i)
                m_pixels[i] = 0.0;
        }

        // Mark pixel as proper to the class instance
        m_linked = 0;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Open FITS image
 *
 * @param fptr FITS file pointer
 *
 * Open FITS image in FITS file. Opening means connecting the FITS file
 * pointer to the image and reading the image and axes dimensions.
 ***************************************************************************/
void GFitsImageFlt::open(__fitsfile* fptr)
{
    // Initialise base class members
    GFitsImage::init_members();

    // Initialise class members
    init_members();

    // Open image
    open_image(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save FITS image
 *
 * Saves the image into the FITS file.
 ***************************************************************************/
void GFitsImageFlt::save(void)
{
    // Save a single precision image
    save_image(__TFLOAT, m_pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close FITS image
 *
 * Closing a FITS image resets the object into its initial state. Closing
 * does NOT save the image into the FITS file. Use the save method for this
 * purpose.
 ***************************************************************************/
void GFitsImageFlt::close(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Free base class members
    GFitsImage::free_members();

    // Initialise base class members
    GFitsImage::init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone FITS image
 *
 * Cloning provides a copy of the FITS file. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GFitsImageFlt* GFitsImageFlt::clone(void) const
{
    // Clone this image
    return new GFitsImageFlt(*this);
}


/*==========================================================================
 =                                                                         =
 =                          GFitsImageFlt friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GFitsImageFlt                 =
 =                                                                         =
 ==========================================================================*/
