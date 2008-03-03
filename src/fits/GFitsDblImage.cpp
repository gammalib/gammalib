/***************************************************************************
 *         GFitsDblImage.cpp  - FITS double precision image class          *
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
#include "GFitsDblImage.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR2 "GFitsDblImage::operator() (int, int)"
#define G_OPEN      "GFitsDblImage::open(fitsfile*)"
#define G_SAVE      "GFitsDblImage::save()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                  GFitsDblImage constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsDblImage::GFitsDblImage() : GFitsImage()
{
    // Initialise class members for clean destruction
    init_members();

    // Set bitpix
    m_bitpix = -64;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param naxis Image dimension (0,1,2,3,4)
 * @param naxes Number of pixels in each dimension
 *
 * Construct instance of GFitsDblImage by specifying the image dimension and
 * the number of pixels in each dimension. Note that this constructor does
 * not allocate any memory for the actual image.
 ***************************************************************************/
GFitsDblImage::GFitsDblImage(int naxis, const int* naxes) : GFitsImage(naxis, naxes)
{
    // Initialise class members for clean destruction
    init_members();

    // Set bitpix
    m_bitpix = -64;

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
 * Construct instance of GFitsDblImage by specifying the image dimension, the
 * number of pixels in each dimension and an image pixel array. If the
 * 'pixels' pointer is valid (i.e. not NULL) the pixel data is indeed copied
 * by this constructor. In case of memory shortage a combination of
 *   GFitsDblImage::GFitsDblImage(int naxis, const int* naxes)
 * and
 *   GFitsDblImage::link(double* pixels)
 * may be prefered to avoid duplication of the pixel data.
 ***************************************************************************/
GFitsDblImage::GFitsDblImage(int naxis, const int* naxes, const double* pixels) 
                                                   : GFitsImage(naxis, naxes)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy pixels into object
    if (m_num_pixels > 0 && pixels != NULL) {
        m_pixels = new double[m_num_pixels];
        memcpy(m_pixels, pixels, m_num_pixels*sizeof(double));
    }

    // Set bitpix
    m_bitpix = -64;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param image FITS image which should be used for construction
 ***************************************************************************/
GFitsDblImage::GFitsDblImage(const GFitsDblImage& image) : GFitsImage(image)
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
GFitsDblImage::~GFitsDblImage()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsDblImage operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param image FITS image to be assigned
 ***************************************************************************/
GFitsDblImage& GFitsDblImage::operator= (const GFitsDblImage& image)
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
double& GFitsDblImage::operator() (const int& ix)
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
const double& GFitsDblImage::operator() (const int& ix) const
{
    // Operator is only valid for 1D images
    if (m_naxis != 1)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 1);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsDblImage*)this)->fetch_pixels();

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
double& GFitsDblImage::operator() (const int& ix, const int& iy)
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
const double& GFitsDblImage::operator() (const int& ix, const int& iy) const
{
    // Operator is only valid for 2D images
    if (m_naxis != 2)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 2);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsDblImage*)this)->fetch_pixels();

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
double& GFitsDblImage::operator() (const int& ix, const int& iy, const int& iz)
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
const double& GFitsDblImage::operator() (const int& ix, const int& iy, const int& iz) const
{
    // Operator is only valid for 3D images
    if (m_naxis != 3)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 3);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsDblImage*)this)->fetch_pixels();

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
double& GFitsDblImage::operator() (const int& ix, const int& iy, const int& iz, const int& it)
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
const double& GFitsDblImage::operator() (const int& ix, const int& iy, const int& iz, const int& it) const
{
    // Operator is only valid for 4D images
    if (m_naxis != 4)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 4);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsDblImage*)this)->fetch_pixels();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2]));

    // Return image pixel
    return m_pixels[offset];
}


/*==========================================================================
 =                                                                         =
 =                       GFitsDblImage public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Open FITS image
 *
 * @param fptr FITS file pointer
 *
 * Open FITS image in FITS file. Opening means connecting the FITS file
 * pointer to the image and reading the image and axes dimensions.
 ***************************************************************************/
void GFitsDblImage::open(__fitsfile* fptr)
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
void GFitsDblImage::link(double* pixels)
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
 * @brief Save FITS image
 *
 * Saves the image into the FITS file.
 ***************************************************************************/
void GFitsDblImage::save(void)
{
    // Save a double precision image
    save_image(__TDOUBLE, m_pixels);

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
void GFitsDblImage::close(void)
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
GFitsDblImage* GFitsDblImage::clone(void) const
{
    // Clone this image
    return new GFitsDblImage(*this);
}


/***********************************************************************//**
 * @brief Set NULL value
 *
 * @param value Value to which NULL pixels should be set
 *
 * NOTE: FOR THE MOMENT THE IMAGE IS NOT RE-LOADED IF THE NULL VALUE IS
 * CHANGED !!!!
 ***************************************************************************/
void GFitsDblImage::set_nullval(const double* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new double;
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


/*==========================================================================
 =                                                                         =
 =                       GFitsDblImage private methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsDblImage::init_members(void)
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
void GFitsDblImage::copy_members(const GFitsDblImage& image)
{
    // Copy attributes
    m_linked = image.m_linked;
    m_anynul = image.m_anynul;

    // Copy pixels
    if (m_linked)
        m_pixels = image.m_pixels;
    else {
        if (m_num_pixels > 0 && image.m_pixels != NULL) {
            m_pixels = new double[m_num_pixels];
            memcpy(m_pixels, image.m_pixels, m_num_pixels*sizeof(double));
        }
    }

    // Copy NULL value
    if (image.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new double;
        *m_nulval = *image.m_nulval;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsDblImage::free_members(void)
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
void GFitsDblImage::fetch_pixels(void)
{
    // Fetch only if there are pixels in image
    if (m_num_pixels > 0) {

        // If we have pixels allocated then free them now
        if (!m_linked) {
            if (m_pixels != NULL) delete [] m_pixels;
        }

        // Allocate fresh image pixels
        m_pixels = new double[m_num_pixels];

        // Case 1: A FITS is attached. Load pixels now.
        if (m_fitsfile.Fptr != NULL)
            load_image(__TDOUBLE, m_pixels, m_nulval, &m_anynul);

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


/*==========================================================================
 =                                                                         =
 =                          GFitsDblImage friends                          =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GFitsDblImage                 =
 =                                                                         =
 ==========================================================================*/
