/***************************************************************************
 *                       GSkymap.hpp - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * This class implements a sky maps. Sky maps are collections of pixels that
 * define a quantity as function of pixel location. Typical quantities are
 * gamma-ray intensities, but may also be the number of measured counts.
 *
 * Sky map pixels may be arranged in a 2-dimensional grid or in a linear
 * 1-dimensional sequence. The link between pixel index and sky direction
 * on the celestial sphere is established using a projection, implemented
 * by the GSkyProjection class. 2-dimensional grids are represented by
 * World Coordinate Systems. All World Coordinate Systems derive from GWcs
 * and are registered in GWcsRegistry. The Healpix pixelisation, implemented
 * by the GHealpix class, is the only 1-dimensional grid that is so far
 * available.
 *
 * Sky map pixels may be accessed by their linear index, by pixel or by
 * sky direction. While the index and pixel access return the sky map value
 * at the pixel centre, the sky direction access operator performs an
 * interpolation to the exact sky direction.
 *
 * Conversion methods exist to convert between the linear index, the pixel
 * and the sky direction:
 *
 *     GSkyPixel pixel = map.inx2pix(index);   // Index to pixel
 *     GSkyDir   dir   = map.inx2dir(index);   // Index to sky direction
 *     GSkyDir   dir   = map.pix2dir(pixel);   // Pixel to sky direction
 *     int       index = map.pix2inx(pixel);   // Pixel to index
 *     int       index = map.dir2inx(dir);     // Sky direction to index
 *     GSkyPixel pixel = map.dir2pix(dir);     // Sky direction to pixel
 *  
 ***************************************************************************/
class GSkymap : public GBase {

	friend GSkymap sqrt(const GSkymap& map);

public:
    // Constructors and destructors
    GSkymap(void);
    explicit GSkymap(const std::string& filename);
    explicit GSkymap(const std::string& coords,
                     const int&         nside,
                     const std::string& order,
                     const int&         nmaps = 1);
    explicit GSkymap(const std::string& wcs,
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
    GSkymap&      operator=(const double& value);
    GSkymap&      operator+=(const GSkymap& map);
    GSkymap&      operator-=(const GSkymap& map);
    GSkymap&      operator*=(const GSkymap& map);
    GSkymap&      operator/=(const GSkymap& map);
    double&       operator()(const int& index, const int& map = 0);
    const double& operator()(const int& index, const int& map = 0) const;
    double&       operator()(const GSkyPixel& pixel, const int& map = 0);
    const double& operator()(const GSkyPixel& pixel, const int& map = 0) const;
    double        operator()(const GSkyDir& dir, const int& map = 0) const;

    // Methods
    void                  clear(void);
    GSkymap*              clone(void) const;
    std::string           classname(void) const;
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
    double                solidangle(const int& index) const;
    double                solidangle(const GSkyPixel& pixel) const;
    bool                  contains(const GSkyDir& dir) const;
    bool                  contains(const GSkyPixel& pixel) const;
    const GSkyProjection* projection(void) const;
    void                  projection(const GSkyProjection& proj);
    const double*         pixels(void) const;
    void                  load(const std::string& filename);
    void                  save(const std::string& filename, bool clobber = false) const;
    void                  read(const GFitsHDU& hdu);
    void                  write(GFits& file) const;
    std::string           print(const GChatter& chatter = NORMAL) const;
    GSkymap        stack_maps(void);

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
    void              read_healpix(const GFitsTable& table);
    void              read_wcs(const GFitsImage& image);
    void              alloc_wcs(const GFitsImage& image);
    GFitsBinTable*    create_healpix_hdu(void) const;
    GFitsImageDouble* create_wcs_hdu(void) const;

    // Private data area
    int             m_num_pixels; //!< Number of pixels (used for pixel allocation)
    int             m_num_maps;   //!< Number of maps (used for pixel allocation)
    int             m_num_x;      //!< Number of pixels in x direction (only 2D)
    int             m_num_y;      //!< Number of pixels in y direction (only 2D)
    GSkyProjection* m_proj;       //!< Pointer to sky projection
    double*         m_pixels;     //!< Pointer to skymap pixels

    // Computation cache
    mutable bool    m_hascache;   //!< Cache is valid
    mutable bool    m_contained;  //!< Direction contained in map
    mutable GSkyDir m_last_dir;   //!< Last sky direction
    mutable int     m_inx1;       //!< Interpolation index 1
    mutable int     m_inx2;       //!< Interpolation index 2
    mutable int     m_inx3;       //!< Interpolation index 3
    mutable int     m_inx4;       //!< Interpolation index 4
    mutable double  m_wgt1;       //!< Interpolation weight 1
    mutable double  m_wgt2;       //!< Interpolation weight 2
    mutable double  m_wgt3;       //!< Interpolation weight 3
    mutable double  m_wgt4;       //!< Interpolation weight 4
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkymap").
 ***************************************************************************/
inline
std::string GSkymap::classname(void) const
{
    return ("GSkymap");
}


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
