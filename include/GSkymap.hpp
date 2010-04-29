/***************************************************************************
 *            GSkymap.hpp  -  Class that implements a sky map              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkymap.hpp
 * @brief GSkymap class definition.
 * @author J. Knodlseder
 */

#ifndef GSKYMAP_HPP
#define GSKYMAP_HPP

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GWcs.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
#include "GFitsHDU.hpp"
#include "GMatrix.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GSkymap
 *
 * @brief GSkymap class interface defintion
 ***************************************************************************/
class GSkymap {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GSkymap& map);

public:
    // Constructors and destructors
    GSkymap(void);
    GSkymap(const std::string& filename);
    GSkymap(const std::string& wcs, const std::string& coords,
            const int& nside, const std::string& order,
            const int nmaps = 1);
    GSkymap(const std::string& wcs, const std::string& coords,
            double const& x, double const& y,
            double const& dx, double const& dy,
            const int& nx, const int& ny, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // Operators
    GSkymap& operator= (const GSkymap& map);

    // 1D pixel methods
    double&       operator() (const int& pixel, const int map = 0);
    const double& operator() (const int& pixel, const int map = 0) const;
    GSkyDir       pix2dir(const int& pix);
    int           dir2pix(GSkyDir dir) const;
    double        omega(const int& pix) const;

    // 2D pixel methods
    double&       operator() (const GSkyPixel& pixel, const int map = 0);
    const double& operator() (const GSkyPixel& pixel, const int map = 0) const;
    GSkyDir       xy2dir(const GSkyPixel& pix);
    GSkyPixel     dir2xy(GSkyDir dir) const;

    // Methods
    void      load(const std::string& filename);
    void      save(const std::string& filename, int clobber = 0);
    void      read(const GFitsHDU* hdu);
    void      write(GFitsHDU* hdu);
    int       npix(void) const;
    int       nx(void) const;
    int       ny(void) const;
    int       nmaps(void) const;

private:
    // Private methods
    void      init_members(void);
    void      alloc_pixels(void);
    void      copy_members(const GSkymap& map);
    void      free_members(void);
    void      set_wcs(const std::string& wcs, const std::string& coords,
                      const double& crval1, const double& crval2,
                      const double& crpix1, const double& crpix2,
                      const double& cdelt1, const double& cdelt2,
                      const GMatrix& cd, const GVector& pv2);
    int       xy2pix(const GSkyPixel& pix) const;
    GSkyPixel pix2xy(const int& pix) const;
    void      read_healpix(const GFitsHDU* hdu);
    GFitsHDU* create_healpix_hdu(void);

    // Private data area
    int     m_num_pixels;   //!< Number of pixels (used for pixel allocation)
    int     m_num_maps;     //!< Number of maps (used for pixel allocation)
    int     m_num_x;        //!< Number of pixels in x direction (only 2D)
    int     m_num_y;        //!< Number of pixels in y direction (only 2D)
    GWcs*   m_wcs;          //!< Pointer to WCS projection
    double* m_pixels;       //!< Pointer to skymap pixels
};

#endif /* GSKYMAP_HPP */
