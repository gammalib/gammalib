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
#include "GWcsHPX.hpp"


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
    GSkymap(const std::string& wcs, const std::string& coordsys, 
            const int& nside, const std::string& ordering,
            const int nmaps = 1);
    GSkymap(const std::string& wcs, const std::string& coordsys, 
            GSkyDir& dir, const int& nlon, const int& nlat,
            const double& dlon, const double& dlat, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // Operators
    GSkymap& operator= (const GSkymap& map);

    // Methods
    void read(const GFitsHDU* hdu);
    void write(const GFits* file);

private:
    // Private methods
    void      init_members(void);
    void      alloc_pixels(void);
    void      copy_members(const GSkymap& map);
    void      free_members(void);
    void      read_healpix(const GFitsHDU* hdu);
    GFitsHDU* create_hdu_healpix(void);

    // Private data area
    int     m_num_pixels; //!< Number of pixels per map
    int     m_num_maps;   //!< Number of maps
    GWcs*   m_wcs;        //!< Pointer to WCS projection
    double* m_pixels;     //!< Pointer to skymap pixels
};

#endif /* GSKYMAP_HPP */
