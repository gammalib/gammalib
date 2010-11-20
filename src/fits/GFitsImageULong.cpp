/***************************************************************************
 *         GFitsImageULong.cpp  - FITS unsigned long image class          *
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
#include "GFitsImageULong.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */
#define G_BITPIX 32                 //!< Defines the number of bits per pixel

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(void) : GFitsImage()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief 1D image constructor
 *
 * @param[in] nx Number of pixels.
 * @param[in] pixels Optional pointer to image pixel array
 *
 * Construct 1D instance by specifying the number of pixels in the image.
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(int nx, const unsigned long* pixels) :
                 GFitsImage(G_BITPIX, nx)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct data
    construct_data(pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief 2D image constructor
 *
 * @param[in] nx Number of pixels in first dimension.
 * @param[in] ny Number of pixels in second dimension.
 * @param[in] pixels Optional pointer to image pixel array
 *
 * Construct 2D image by specifying the number of pixels in each dimension.
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(int nx, int ny, const unsigned long* pixels) :
                 GFitsImage(G_BITPIX, nx, ny)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct data
    construct_data(pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief 3D image constructor
 *
 * @param[in] nx Number of pixels in first dimension.
 * @param[in] ny Number of pixels in second dimension.
 * @param[in] nz Number of pixels in third dimension.
 * @param[in] pixels Optional pointer to image pixel array
 *
 * Construct 3D image by specifying the number of pixels in each dimension.
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(int nx, int ny, int nz,
                                 const unsigned long* pixels) :
                 GFitsImage(G_BITPIX, nx, ny, nz)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct data
    construct_data(pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief 4D image constructor
 *
 * @param[in] nx Number of pixels in first dimension.
 * @param[in] ny Number of pixels in second dimension.
 * @param[in] nz Number of pixels in third dimension.
 * @param[in] nt Number of pixels in forth dimension.
 * @param[in] pixels Optional pointer to image pixel array
 *
 * Construct 4D image by specifying the number of pixels in each dimension.
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(int nx, int ny, int nz, int nt,
                                 const unsigned long* pixels) :
                 GFitsImage(G_BITPIX, nx, ny, nz, nt)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct data
    construct_data(pixels);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Pixel array constructor
 *
 * @param[in] naxis Image dimension (0,1,2,3,4).
 * @param[in] naxes Number of pixels in each dimension.
 * @param[in] pixels Optional pointer to image pixel array
 *
 * Construct instance of GFitsImageULong by specifying the image dimension and
 * the number of pixels in each dimension. Note that this constructor does
 * not allocate any memory for the actual image.
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(int naxis, const int* naxes,
                                 const unsigned long* pixels) :
                 GFitsImage(G_BITPIX, naxis, naxes)
{
    // Initialise class members for clean destruction
    init_members();

    // If there are pixels then allocate array
    if (m_num_pixels > 0) {

        // Allocate data
        alloc_data();

        // If no pixel array has been specified then simply initialise data
        if (pixels == NULL)
            init_data();

        // ... otherwise copy pixels
        else {
            for (int i = 0; i < m_num_pixels; ++i)
                m_pixels[i] = pixels[i];
        }

    } // endif: there are pixels in image

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] image FITS image.
 ***************************************************************************/
GFitsImageULong::GFitsImageULong(const GFitsImageULong& image) :
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
GFitsImageULong::~GFitsImageULong(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] image FITS image.
 ***************************************************************************/
GFitsImageULong& GFitsImageULong::operator= (const GFitsImageULong& image)
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
 * @brief Image pixel access operator
 *
 * @param[in] ix Pixel index (starting from 0).
 *
 * Provides access to an image pixel. No range checking is performed.
 * Use the at(ix) method if range checking is required.
 ***************************************************************************/
unsigned long& GFitsImageULong::operator() (const int& ix)
{
    // If image pixels are not available then fetch them now
    if (m_pixels == NULL) fetch_data();

    // Return image pixel
    return m_pixels[ix];
}


/***********************************************************************//**
 * @brief 2D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 *
 * Provides access to a pixel of a 2D image. No range checking or image
 * dimension verification is performed. Use the at(ix,iy) method if range
 * checking and image dimension verification is required.
 ***************************************************************************/
unsigned long& GFitsImageULong::operator() (const int& ix, const int& iy)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Calculate pixel offset
    int offset = ix + iy * m_naxes[0];

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 3D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 *
 * Provides access to a pixel of a 3D image. No range checking or image
 * dimension verification is performed. Use the at(ix,iy) method if range
 * checking and image dimension verification is required.
 ***************************************************************************/
unsigned long& GFitsImageULong::operator() (const int& ix, const int& iy,
                                            const int& iz)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + iz * m_naxes[1]);

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 4D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 * @param[in] it Pixel index in forth dimension (starting from 0).
 *
 * Provides access to a pixel of a 4D image. No range checking or image
 * dimension verification is performed. Use the at(ix,iy) method if range
 * checking and image dimension verification is required.
 ***************************************************************************/
unsigned long& GFitsImageULong::operator() (const int& ix, const int& iy,
                                            const int& iz, const int& it)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2]));

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief Image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index (starting from 0).
 *
 * Provides access to an image pixel. No range checking is performed.
 * Use the at(ix) method if range checking is required.
 ***************************************************************************/
const unsigned long& GFitsImageULong::operator() (const int& ix) const
{
    // If image pixels are not available then fetch them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Return image pixel
    return m_pixels[ix];
}


/***********************************************************************//**
 * @brief 2D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 *
 * Provides access to a pixel of a 2D image. No range checking or image
 * dimension verification is performed. Use the at(ix,iy) method if range
 * checking and image dimension verification is required.
 ***************************************************************************/
const unsigned long& GFitsImageULong::operator() (const int& ix, const int& iy) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Calculate pixel offset
    int offset = ix + iy * m_naxes[0];

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 3D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 *
 * Provides access to a pixel of a 3D image. No range checking or image
 * dimension verification is performed. Use the at(ix,iy) method if range
 * checking and image dimension verification is required.
 ***************************************************************************/
const unsigned long& GFitsImageULong::operator() (const int& ix, const int& iy,
                                                  const int& iz) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + iz * m_naxes[1]);

    // Return image pixel
    return m_pixels[offset];
}


/***********************************************************************//**
 * @brief 4D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 * @param[in] it Pixel index in forth dimension (starting from 0).
 *
 * Provides access to a pixel of a 4D image. No range checking or image
 * dimension verification is performed. Use the at(ix,iy) method if range
 * checking and image dimension verification is required.
 ***************************************************************************/
const unsigned long& GFitsImageULong::operator() (const int& ix, const int& iy,
                                                  const int& iz, const int& it) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Calculate pixel offset
    int offset = ix + m_naxes[0] * (iy + m_naxes[1] * (iz + it *  m_naxes[2]));

    // Return image pixel
    return m_pixels[offset];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Image pixel access operator
 *
 * @param[in] ix Pixel index (starting from 0).
 *
 * Provides access to a pixel of an image including range checking. Note that
 * this method does not necessarily restrict do 1D images but generally
 * applies to all image dimensions.
 ***************************************************************************/
unsigned long& GFitsImageULong::at(const int& ix)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix)]);
}


/***********************************************************************//**
 * @brief 2D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
unsigned long& GFitsImageULong::at(const int& ix, const int& iy)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix,iy)]);
}


/***********************************************************************//**
 * @brief 3D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
unsigned long& GFitsImageULong::at(const int& ix, const int& iy, const int& iz)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix,iy,iz)]);
}


/***********************************************************************//**
 * @brief 4D image pixel access operator
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 * @param[in] it Pixel index in forth dimension (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
unsigned long& GFitsImageULong::at(const int& ix, const int& iy, const int& iz,
                                   const int& it)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix,iy,iz,it)]);
}


/***********************************************************************//**
 * @brief Image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
const unsigned long& GFitsImageULong::at(const int& ix) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix)]);
}


/***********************************************************************//**
 * @brief 2D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
const unsigned long& GFitsImageULong::at(const int& ix, const int& iy) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix,iy)]);
}


/***********************************************************************//**
 * @brief 3D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
const unsigned long& GFitsImageULong::at(const int& ix, const int& iy,
                                         const int& iz) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix,iy,iz)]);
}


/***********************************************************************//**
 * @brief 4D image pixel access operator (const variant)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 * @param[in] it Pixel index in forth dimension (starting from 0).
 *
 * Provides access to a pixel of an image including range checking.
 ***************************************************************************/
const unsigned long& GFitsImageULong::at(const int& ix, const int& iy,
                                         const int& iz, const int& it) const
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) ((GFitsImageULong*)this)->fetch_data();

    // Return image pixel
    return (m_pixels[offset(ix,iy,iz,it)]);
}


/***********************************************************************//**
 * @brief Return value of image pixel
 *
 * @param[in] ix Pixel index (starting from 0).
 *
 * Returns the value of an image pixel as double precision floating point
 * value. This method performs range checking.
 ***************************************************************************/
double GFitsImageULong::pixel(const int& ix) const
{
    // Return pixel value
    return (double(this->at(ix)));
}


/***********************************************************************//**
 * @brief Return value of image pixel (2D version)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 *
 * Returns the value of an image pixel as double precision floating point
 * value. This method performs range checking.
 ***************************************************************************/
double GFitsImageULong::pixel(const int& ix, const int& iy) const
{
    // Return pixel value
    return (double(this->at(ix,iy)));
}


/***********************************************************************//**
 * @brief Return value of image pixel (3D version)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 *
 * Returns the value of an image pixel as double precision floating point
 * value. This method performs range checking.
 ***************************************************************************/
double GFitsImageULong::pixel(const int& ix, const int& iy, const int& iz) const
{
    // Return pixel value
    return (double(this->at(ix,iy,iz)));
}


/***********************************************************************//**
 * @brief Return value of image pixel (4D version)
 *
 * @param[in] ix Pixel index in first dimension (starting from 0).
 * @param[in] iy Pixel index in second dimension (starting from 0).
 * @param[in] iz Pixel index in third dimension (starting from 0).
 * @param[in] it Pixel index in forth dimension (starting from 0).
 *
 * Returns the value of an image pixel as double precision floating point
 * value. This method performs range checking.
 ***************************************************************************/
double GFitsImageULong::pixel(const int& ix, const int& iy, const int& iz,
                              const int& it) const
{
    // Return pixel value
    return (double(this->at(ix,iy,iz,it)));
}


/***********************************************************************//**
 * @brief Return pointer to image pixel
 ***************************************************************************/
void* GFitsImageULong::pixels(void)
{
    // If image pixels are not available then allocate them now
    if (m_pixels == NULL) fetch_data();

    // Return
    return m_pixels;
}


/***********************************************************************//**
 * @brief Return image type
 ***************************************************************************/
int GFitsImageULong::type(void) const
{
    // Return type
    return __TULONG;
}


/***********************************************************************//**
 * @brief Clone FITS image
 *
 * Cloning provides a copy of the FITS file. Cloning is used to allocate
 * derived classes into a base class pointer.
 ***************************************************************************/
GFitsImageULong* GFitsImageULong::clone(void) const
{
    // Clone this image
    return new GFitsImageULong(*this);
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsImageULong::init_members(void)
{
    // Initialise members
    m_bitpix = G_BITPIX;
    m_pixels = NULL;
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
void GFitsImageULong::copy_members(const GFitsImageULong& image)
{
    // Fetch column data if not yet fetched. The casting circumvents the
    // const correctness
    bool not_loaded = (image.m_pixels == NULL);
    if (not_loaded) ((GFitsImageULong*)(&image))->fetch_data();

    // Copy pixels
    if (m_num_pixels > 0 && image.m_pixels != NULL) {
        m_pixels = new unsigned long[m_num_pixels];
        for (int i = 0; i < m_num_pixels; ++i)
            m_pixels[i] = image.m_pixels[i];
    }

    // Copy NULL value
    alloc_nulval(image.m_nulval);

    // Small memory option: release column if it was fetch above
    #if defined(G_SMALL_MEMORY)
    if (not_loaded) ((GFitsImageULong*)(&image))->release_data();
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsImageULong::free_members(void)
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
 * @brief Allocate data
 ***************************************************************************/
void GFitsImageULong::alloc_data(void)
{
    // Release any existing data
    release_data();

    // Allocate new data
    if (m_num_pixels > 0)
        m_pixels = new unsigned long[m_num_pixels];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise data
 ***************************************************************************/
void GFitsImageULong::init_data(void)
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
void GFitsImageULong::release_data(void)
{
    // Free any existing memory
    if (m_pixels != NULL) delete [] m_pixels;

    // Mark pointer as free
    m_pixels = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct data from array
 *
 * @param[in] pixels Optional pointer to pixel array.
 *
 * Initialise pixel data from an optional pixel array. If the pointer is
 * NULL the image is simply initialised. This method supports all
 * constructors.
 ***************************************************************************/
void GFitsImageULong::construct_data(const unsigned long* pixels)
{
    // If there are pixels then allocate array
    if (m_num_pixels > 0) {

        // Allocate data
        alloc_data();

        // If no pixel array has been specified then simply initialise data
        if (pixels == NULL)
            init_data();

        // ... otherwise copy pixels
        else {
            for (int i = 0; i < m_num_pixels; ++i)
                m_pixels[i] = pixels[i];
        }

    } // endif: there are pixels in image

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocates nul value
 ***************************************************************************/
void GFitsImageULong::alloc_nulval(const void* value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // If we have valid value, allocate and set nul value
    if (value != NULL) {
        m_nulval  = new unsigned long;
        *m_nulval = *((unsigned long*)value);
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
