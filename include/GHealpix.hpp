/***************************************************************************
 *             GHealpix.hpp  -  Healpix sky representation class           *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GHealpix.hpp
 * @brief GHealpix class definition.
 * @author J. Knodlseder
 */

#ifndef GHEALPIX_HPP
#define GHEALPIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHDU.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief GHealpix class interface defintion
 ***************************************************************************/
class GHealpix {

public:
    // Constructors and destructors
    GHealpix();
    GHealpix(const GHealpix& pixels);
    virtual ~GHealpix();

    // Operators
    GHealpix& operator= (const GHealpix& pixels);

    // Methods
    void   load(const GFitsHDU* hdu);
    int    nside(void) const;
    int    num_pixels(void) const;
    double omega(void) const;
    
private:
    // Private methods
    void      init_members(void);
    void      copy_members(const GHealpix& pixels);
    void      free_members(void);
    GHealpix* clone(void) const;

    // Private data area
    int     m_nside;        //!< Number of divisions of each base pixel (1,2,..)
    int     m_order;        //!< Ordering (0=ring, 1=nested)
    int     m_coordsys;     //!< Coordinate system (0=equatorial, 1=galactic)
    int     m_num_pixels;   //!< Number of pixels
    int     m_size_pixels;  //!< Vector size of each pixel
    double* m_pixels;       //!< Pixel array
    double* m_lon;          //!< RA or GLON values for each pixel
    double* m_lat;          //!< DEC or GLAT values for each pixel
};

#endif /* GHEALPIX_HPP */
