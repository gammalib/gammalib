/***************************************************************************
 *            GSkymap.hpp  -  Class that implements a sky map              *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2010 by Jurgen Knodlseder                   *
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
#include "GSkyDir.hpp"
#include "GWcs.hpp"
#include "GFitsHDU.hpp"


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
    GSkymap(const std::string& wcs, const std::string& coords,
            const int& nside, const std::string& order,
            const int nmaps = 1);
    GSkymap(const std::string& wcs, const std::string& coords,
            GSkyDir& dir, const int& nlon, const int& nlat,
            const double& dlon, const double& dlat, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // Pixel access operators
    double&       operator() (int pixel, int element = 0);
    const double& operator() (int pixel, int element = 0) const;

    // Operators
    GSkymap& operator= (const GSkymap& map);

    // Methods
    void    load(const std::string& filename);
    void    save(const std::string& filename, int clobber = 0);
    void    read(const GFitsHDU* hdu);
    void    write(GFitsHDU* hdu);
    GSkyDir pix2dir(const int& ipix);
    int     dir2pix(GSkyDir dir) const;
    double  omega(const int& pix) const;
    int     npix(void) const;

private:
    // Private methods
    void      init_members(void);
    void      alloc_pixels(void);
    void      copy_members(const GSkymap& map);
    void      free_members(void);
    void      read_healpix(const GFitsHDU* hdu);
    GFitsHDU* create_healpix_hdu(void);

    // Private data area
    int     m_num_pixels;   //!< Number of pixels
    int     m_num_maps;     //!< Number of maps
    GWcs*   m_wcs;          //!< Pointer to WCS projection
    double* m_pixels;       //!< Pointer to skymap pixels
};

#endif /* GSKYMAP_HPP */
