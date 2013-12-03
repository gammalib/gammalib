/***************************************************************************
 *                       GSkymap.hpp - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GSkymap.hpp
 * @brief Sky map class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSKYMAP_HPP
#define GSKYMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
#include "GSkyProjection.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageDouble.hpp"
#include "GMatrix.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GSkymap
 *
 * @brief Sky map class
 *
 * This class implements a set of sky maps. Skymaps may either be present in
 * the HEALPix format or in any (supported) WCS format. Skymap pixels are
 * stored in a double precision array indexed by (x,y,map), with the x
 * axis being the most rapidely varying axis. 
 * 
 * Skymap pixels may be either accessed via index operators (i,map), where i 
 * runs over the (x,y) direction of the map, or via pixel operators (pixel,map),
 * where pixel implements a 2D (x,y) pixel direction. The first operator is
 * the preferred access method for HEALPix maps while the second operator is
 * the preferred access method for WCS maps. Conversion methods between the
 * index or sky pixel and the true physical sky direction are provided by the
 * pix2dir() and dir2pix() methods.
 *
 * @todo Rewrite class description.
 ***************************************************************************/
class GSkymap : public GBase {

public:
    // Constructors and destructors
    GSkymap(void);
    explicit GSkymap(const std::string& filename);
    explicit GSkymap(const std::string& coords,
                     const int&         nside,
                     const std::string& order,
                     const int&         nmaps = 1);
    explicit GSkymap(const std::string& proj,
                     const std::string& coords,
                     const double&      x,
                     const double&      y,
                     const double&      dx,
                     const double&      dy,
                     const int&         nx,
                     const int&         ny,
                     const int&         nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // Operators
    GSkymap&      operator=(const GSkymap& map);
    double&       operator()(const int& index, const int& map = 0);
    const double& operator()(const int& index, const int& map = 0) const;
    double&       operator()(const GSkyPixel& pixel, const int& map = 0);
    const double& operator()(const GSkyPixel& pixel, const int& map = 0) const;
    double        operator()(const GSkyDir& dir, const int& map = 0) const;

    // Methods
    void                  clear(void);
    GSkymap*              clone(void) const;
    const int&            npix(void) const;
    const int&            nx(void) const;
    const int&            ny(void) const;
    const int&            nmaps(void) const;
    GSkyPixel             inx2pix(const int& index) const;
    GSkyDir               inx2dir(const int& index) const;
    GSkyDir               pix2dir(const GSkyPixel& pixel) const;
    int                   pix2inx(const GSkyPixel& pixel) const;
    int                   dir2inx(const GSkyDir& dir) const;
    GSkyPixel             dir2pix(const GSkyDir& dir) const;
    double                omega(const int& index) const;
    double                omega(const GSkyPixel& pixel) const;
    bool                  contains(const GSkyDir& dir) const;
    bool                  contains(const GSkyPixel& pixel) const;
    const GSkyProjection* projection(void) const;
    void                  projection(const GSkyProjection& proj);
    const double*         pixels(void) const;
    void                  load(const std::string& filename);
    void                  save(const std::string& filename, bool clobber = false) const;
    void                  read(const GFitsHDU* hdu);
    void                  write(GFits* file) const;
    std::string           print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void              init_members(void);
    void              alloc_pixels(void);
    void              copy_members(const GSkymap& map);
    void              free_members(void);
    void              set_wcs(const std::string& wcs, const std::string& coords,
                              const double& crval1, const double& crval2,
                              const double& crpix1, const double& crpix2,
                              const double& cdelt1, const double& cdelt2,
                              const GMatrix& cd, const GVector& pv2);
    void              read_healpix(const GFitsTable* hdu);
    void              read_wcs(const GFitsImage* hdu);
    void              alloc_wcs(const GFitsImage* hdu);
    GFitsBinTable*    create_healpix_hdu(void) const;
    GFitsImageDouble* create_wcs_hdu(void) const;

    // Private data area
    int             m_num_pixels; //!< Number of pixels (used for pixel allocation)
    int             m_num_maps;   //!< Number of maps (used for pixel allocation)
    int             m_num_x;      //!< Number of pixels in x direction (only 2D)
    int             m_num_y;      //!< Number of pixels in y direction (only 2D)
    GSkyProjection* m_proj;       //!< Pointer to sky projection
    double*         m_pixels;     //!< Pointer to skymap pixels
};


/***********************************************************************//**
 * @brief Returns number of pixels
 *
 * @return Number of pixels in one sky map.
 *
 * Returns the number of pixels in one sky map.
 ***************************************************************************/
inline
const int& GSkymap::npix(void) const
{
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns number of pixels in x coordinate
 *
 * @return Number of pixels in the X direction.
 *
 * Returns the number of pixels in the X direction. If the sky map is a
 * one dimensional array (which is the case for the Healpix projection),
 * the method returns 0.
 ***************************************************************************/
inline
const int& GSkymap::nx(void) const
{
    return m_num_x;
}


/***********************************************************************//**
 * @brief Returns number of pixels in y coordinate
 *
 * @return Number of pixels in the Y direction.
 *
 * Returns the number of pixels in the Y direction. If the sky map is a
 * one dimensional array (which is the case for the Healpix projection),
 * the method returns 0.
 ***************************************************************************/
inline
const int& GSkymap::ny(void) const
{
    return m_num_y;
}


/***********************************************************************//**
 * @brief Returns number of maps
 *
 * @return Number of maps in the sky map object.
 ***************************************************************************/
inline
const int& GSkymap::nmaps(void) const
{
    return m_num_maps;
}


/***********************************************************************//**
 * @brief Returns pointer to sky projection
 *
 * @return Pointer to sky projection (NULL if no projection is defined).
 ***************************************************************************/
inline
const GSkyProjection* GSkymap::projection(void) const
{
    return m_proj;
}


/***********************************************************************//**
 * @brief Returns pointer to pixel data
 *
 * @return Pointer to pixel data.
 ***************************************************************************/
inline
const double* GSkymap::pixels(void) const
{
    return m_pixels;
}

#endif /* GSKYMAP_HPP */
