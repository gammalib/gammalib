/***************************************************************************
 *               GFitsImage.hpp - Abstract FITS image base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsImage.hpp
 * @brief Abstract FITS image base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSIMAGE_HPP
#define GFITSIMAGE_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHDU.hpp"


/***********************************************************************//**
 * @class GFitsImage
 *
 * @brief Abstract FITS image base class
 *
 * This class defines the abstract interface for a FITS image.
 ***************************************************************************/
class GFitsImage : public GFitsHDU {

public:
    // Constructors and destructors
    GFitsImage(void);
    GFitsImage(const int& bitpix, const int& nx);
    GFitsImage(const int& bitpix, const int& nx, const int& ny);
    GFitsImage(const int& bitpix, const int& nx, const int& ny, const int& nz);
    GFitsImage(const int& bitpix, const int& nx, const int& ny, const int& nz, const int& nt);
    GFitsImage(const int& bitpix, const int& naxis, const int* naxes);
    GFitsImage(const GFitsImage& image);
    virtual ~GFitsImage(void);

    // Operators
    GFitsImage& operator=(const GFitsImage& image);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GFitsImage* clone(void) const = 0;
    virtual void*       pixels(void) = 0;
    virtual double      pixel(const int& ix) const = 0;
    virtual double      pixel(const int& ix, const int& iy) const = 0;
    virtual double      pixel(const int& ix, const int& iy, const int& iz) const = 0;
    virtual double      pixel(const int& ix, const int& iy, const int& iz, const int& it) const = 0;
    virtual int         type(void) const = 0;

    // Implemented pure virtual methods
    HDUType exttype(void) const;

    // Base class methods
    const int&  size(void) const;
    const int&  bitpix(void) const;
    const int&  naxis(void) const;
    int         naxes(const int& axis) const;
    const int&  anynul(void) const;
    void        nulval(const void* value);
    const void* nulval(void) const;
    std::string print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GFitsImage& image);
    void  free_members(void);
    void  data_open(void* vptr);
    void  data_save(void);
    void  data_close(void);
    void  data_connect(void* vptr);
    void  init_image_header(void);
    void  open_image(void* vptr);
    void  load_image(int datatype, const void* pixels,
                     const void* nulval, int* anynul);
    void  save_image(int datatype, const void* pixels);
    void  fetch_data(void);
    int   offset(const int& ix) const;
    int   offset(const int& ix, const int& iy) const;
    int   offset(const int& ix, const int& iy, const int& iz) const;
    int   offset(const int& ix, const int& iy, const int& iz, const int& it) const;

    // Pure virtual protected methods
    virtual void  alloc_data(void) = 0;
    virtual void  init_data(void) = 0;
    virtual void  release_data(void) = 0;
    virtual void  alloc_nulval(const void* value) = 0;
    virtual void* ptr_data(void) = 0;
    virtual void* ptr_nulval(void) = 0;

    // Protected data area
    int   m_bitpix;      //!< Number of Bits/pixel
    int   m_naxis;       //!< Image dimension
    long* m_naxes;       //!< Number of pixels in each dimension
    int   m_num_pixels;  //!< Number of image pixels
    int   m_anynul;      //!< Number of NULLs encountered
};


/***********************************************************************//**
 * @brief Return extension type
 *
 * @return Extension type (HT_IMAGE).
 ***************************************************************************/
inline
GFitsHDU::HDUType GFitsImage::exttype(void) const
{
    return (HT_IMAGE);
}


/***********************************************************************//**
 * @brief Return size of pixel array
 ***************************************************************************/
inline
const int& GFitsImage::size(void) const
{
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Return number of Bits per pixel (negative=floating point)
 ***************************************************************************/
inline
const int& GFitsImage::bitpix(void) const
{
    return m_bitpix;
}


/***********************************************************************//**
 * @brief Return dimension of image
 ***************************************************************************/
inline
const int& GFitsImage::naxis(void) const
{
    return m_naxis;
}


/***********************************************************************//**
 * @brief Return number of nul values envountered during loading
 ***************************************************************************/
inline
const int& GFitsImage::anynul(void) const
{
    return m_anynul;
}


/***********************************************************************//**
 * @brief Return nul value
 ***************************************************************************/
inline
const void* GFitsImage::nulval(void) const
{
    return (const_cast<GFitsImage*>(this)->ptr_nulval());
}

#endif /* GFITSIMAGE_HPP */
