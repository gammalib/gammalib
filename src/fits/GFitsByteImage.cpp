/***************************************************************************
 *              GFitsByteImage.cpp  - FITS Byte image class                *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2010 by Jurgen Knodlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "GException.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsByteImage.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OPERATOR2                    "GFitsByteImage::operator() (int,int)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                          Constructors/destructors                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsByteImage::GFitsByteImage(void) : GFitsImage()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param naxis Image dimension (0,1,2,3,4)
 * @param naxes Number of pixels in each dimension
 *
 * Construct instance of GFitsByteImage by specifying the image dimension and
 * the number of pixels in each dimension. Note that this constructor does
 * not allocate any memory for the actual image.
 ***************************************************************************/
GFitsByteImage::GFitsByteImage(int naxis, const int* naxes) :
                GFitsImage(8, naxis, naxes)
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate and initialise pixels
    if (m_num_pixels > 0) {
        m_pixels = new unsigned char[m_num_pixels];
        for (int i = 0; i < m_num_pixels; ++i)
            m_pixels[i] = 0;
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param image FITS image which should be used for construction
 ***************************************************************************/
GFitsByteImage::GFitsByteImage(const GFitsByteImage& image) :
                GFitsImage(image)
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
GFitsByteImage::~GFitsByteImage(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param image FITS image to be assigned
 ***************************************************************************/
GFitsByteImage& GFitsByteImage::operator= (const GFitsByteImage& image)
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
unsigned char& GFitsByteImage::operator() (const int& ix)
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
const unsigned char& GFitsByteImage::operator() (const int& ix) const
{
    // Operator is only valid for 1D images
    //if (m_naxis != 1)
    //    throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 1);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsByteImage*)this)->fetch_pixels();

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
unsigned char& GFitsByteImage::operator() (const int& ix, const int& iy)
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
const unsigned char& GFitsByteImage::operator() (const int& ix, const int& iy) const
{
    // Operator is only valid for 2D images
    if (m_naxis != 2)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 2);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsByteImage*)this)->fetch_pixels();

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
unsigned char& GFitsByteImage::operator() (const int& ix, const int& iy, const int& iz)
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
const unsigned char& GFitsByteImage::operator() (const int& ix, const int& iy, const int& iz) const
{
    // Operator is only valid for 3D images
    if (m_naxis != 3)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 3);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsByteImage*)this)->fetch_pixels();

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
unsigned char& GFitsByteImage::operator() (const int& ix, const int& iy, const int& iz, const int& it)
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
const unsigned char& GFitsByteImage::operator() (const int& ix, const int& iy, const int& iz, const int& it) const
{
    // Operator is only valid for 4D images
    if (m_naxis != 4)
        throw GException::fits_wrong_image_operator(G_OPERATOR2, m_naxis, 4);

    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsByteImage*)this)->fetch_pixels();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2]));

    // Return image pixel
    return m_pixels[offset];
}


/*==========================================================================
 =                                                                         =
 =                       GFitsByteImage public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set NULL value
 *
 * @param value Value to which NULL pixels should be set
 *
 * NOTE: FOR THE MOMENT THE IMAGE IS NOT RE-LOADED IF THE NULL VALUE IS
 * CHANGED !!!!
 ***************************************************************************/
void GFitsByteImage::nulval(const unsigned char* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new unsigned char;
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
void* GFitsByteImage::pixels(void)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_pixels();

    // Return
    return m_pixels;
}


/***********************************************************************//**
 * @brief Clone FITS image
 *
 * Cloning provides a copy of the FITS file. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GFitsByteImage* GFitsByteImage::clone(void) const
{
    // Clone this image
    return new GFitsByteImage(*this);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsByteImage::init_members(void)
{
    // Initialise members
    m_pixels = NULL;
    m_nulval = NULL;
    m_anynul = 0;

    // Initialise image
    //init_image();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param image FITS image to copy
 *
 * Copy class members.
 ***************************************************************************/
void GFitsByteImage::copy_members(const GFitsByteImage& image)
{
    // Copy attributes
    m_anynul = image.m_anynul;

    // Copy pixels
    if (m_num_pixels > 0) {
        m_pixels = new unsigned char[m_num_pixels];
        for (int i = 0; i < m_num_pixels; ++i)
            m_pixels[i] = image.m_pixels[i];
    }

    // Copy NULL value
    if (image.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new unsigned char;
        *m_nulval = *image.m_nulval;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsByteImage::free_members(void)
{
    // Free memory
    if (m_pixels != NULL) delete [] m_pixels;
    if (m_nulval != NULL) delete m_nulval;

    // Mark memory as free
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
void GFitsByteImage::fetch_pixels(void)
{
    // Fetch only if there are pixels in image
    if (m_num_pixels > 0) {

        // If we have pixels allocated then free them now
        if (m_pixels != NULL) delete [] m_pixels;

        // Allocate fresh image pixels
        m_pixels = new unsigned char[m_num_pixels];

        // Case 1: A FITS is attached. Load pixels now.
        if (FPTR(m_fitsfile)->Fptr != NULL)
            load_image(__TBYTE, m_pixels, m_nulval, &m_anynul);

        // Case 2: No FITS file is attached. Allocate pixels and set them all to 0.0
        else {
            for (int i = 0; i < m_num_pixels; ++i)
                m_pixels[i] = 0.0;
        }

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate data
 ***************************************************************************/
void GFitsByteImage::alloc_data(void)
{
    // Release any existing data
    release_data();

    // Allocate new data
    if (m_num_pixels > 0)
        m_pixels = new unsigned char[m_num_pixels];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise data
 ***************************************************************************/
void GFitsByteImage::init_data(void)
{
    // Initialise data if they exist
    if (m_pixels != NULL) {
        for (int i = 0; i < m_num_pixels; ++i)
            m_pixels[i] = 0;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Release data
 ***************************************************************************/
void GFitsByteImage::release_data(void)
{
    // Free any existing memory
    if (m_pixels != NULL) delete [] m_pixels;

    // Mark pointer as free
    m_pixels = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocates nul value
 ***************************************************************************/
void GFitsByteImage::alloc_nulval(const void* value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // If we have valid value, allocate and set nul value
    if (value != NULL) {
        m_nulval  = new unsigned char;
        *m_nulval = *((char*)value);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return image type
 ***************************************************************************/
int GFitsByteImage::type(void) const
{
    // Return type
    return __TBYTE;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
