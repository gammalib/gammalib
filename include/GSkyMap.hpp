/***************************************************************************
 *                       GSkyMap.hpp - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GSkyMap.hpp
 * @brief Sky map class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSKYMAP_HPP
#define GSKYMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
#include "GSkyProjection.hpp"
#include "GBilinear.hpp"
#include "GNdarray.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFits;
class GFitsHDU;
class GFitsTable;
class GFitsBinTable;
class GFitsImage;
class GFitsImageDouble;
class GMatrix;
class GVector;
class GSkyRegion;
class GSkyRegionCircle;
class GSkyRegions;


/***********************************************************************//**
 * @class GSkyMap
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
class GSkyMap : public GBase {

	friend GSkyMap sqrt(const GSkyMap& map);
	friend GSkyMap log(const GSkyMap& map);
	friend GSkyMap log10(const GSkyMap& map);
	friend GSkyMap abs(const GSkyMap& map);
	friend GSkyMap sign(const GSkyMap& map);
    friend GSkyMap clip(const GSkyMap& map, const double& thresh);

public:
    // Constructors and destructors
    GSkyMap(void);
    explicit GSkyMap(const GFilename& filename);
    explicit GSkyMap(const GFitsHDU& hdu);
    GSkyMap(const std::string& coords,
            const int&         nside,
            const std::string& order,
            const int&         nmaps = 1);
    GSkyMap(const std::string& wcs,
            const std::string& coords,
            const double&      x,
            const double&      y,
            const double&      dx,
            const double&      dy,
            const int&         nx,
            const int&         ny,
            const int&         nmaps = 1);
    GSkyMap(const GSkyMap& map);
    virtual ~GSkyMap(void);

    // Operators
    GSkyMap&      operator=(const GSkyMap& map);
    GSkyMap&      operator=(const double& value);
    GSkyMap&      operator+=(const GSkyMap& map);
    GSkyMap&      operator+=(const double& value);
    GSkyMap&      operator-=(const GSkyMap& map);
    GSkyMap&      operator-=(const double& value);
    GSkyMap&      operator*=(const GSkyMap& map);
    GSkyMap&      operator*=(const double& factor);
    GSkyMap&      operator/=(const GSkyMap& map);
    GSkyMap&      operator/=(const double& factor);
    GSkyMap       operator+(const GSkyMap& map) const;
    GSkyMap       operator-(const GSkyMap& map) const;
    GSkyMap       operator*(const GSkyMap& map) const;
    GSkyMap       operator/(const GSkyMap& map) const;
    double&       operator()(const int& index, const int& map = 0);
    const double& operator()(const int& index, const int& map = 0) const;
    double&       operator()(const GSkyPixel& pixel, const int& map = 0);
    const double& operator()(const GSkyPixel& pixel, const int& map = 0) const;
    double        operator()(const GSkyDir& dir, const int& map = 0) const;

    // Methods
    void                    clear(void);
    GSkyMap*                clone(void) const;
    std::string             classname(void) const;
    bool                    is_empty(void) const;
    const int&              npix(void) const;
    const int&              nx(void) const;
    const int&              ny(void) const;
    const int&              nmaps(void) const;
    void                    nmaps(const int& nmaps);
    const std::vector<int>& shape(void) const;
    void                    shape(const int& s1);
    void                    shape(const int& s1, const int& s2);
    void                    shape(const int& s1, const int& s2, const int& s3);
    void                    shape(const std::vector<int>& shape);
    int                     ndim(void) const;
    GSkyPixel               inx2pix(const int& index) const;
    GSkyDir                 inx2dir(const int& index) const;
    GSkyDir                 pix2dir(const GSkyPixel& pixel) const;
    int                     pix2inx(const GSkyPixel& pixel) const;
    int                     dir2inx(const GSkyDir& dir) const;
    GSkyPixel               dir2pix(const GSkyDir& dir) const;
    GNdarray                counts(void) const;
    double                  flux(const int& index, const int& map = 0) const;
    double                  flux(const GSkyPixel& pixel, const int& map = 0) const;
    GNdarray                flux(void) const;
    double                  solidangle(const int& index) const;
    double                  solidangle(const GSkyPixel& pixel) const;
    bool                    contains(const GSkyDir& dir) const;
    bool                    contains(const GSkyPixel& pixel) const;
    bool                    overlaps(const GSkyRegion& region) const;
    void                    smooth(const std::string& kernel, const double& par);
    const GSkyProjection*   projection(void) const;
    void                    projection(const GSkyProjection& proj);
    const double*           pixels(void) const;
    GSkyMap                 extract(const int& map, const int& nmaps = 1) const;
    GSkyMap                 extract(const int& startx, const int& stopx,
                                    const int& starty, const int& stopy) const;
    GSkyMap                 extract(const GSkyRegions& inclusions) const;
    void                    stack_maps(void);
    void                    load(const GFilename& filename);
    void                    save(const GFilename& filename,
                                 const bool&      clobber = false) const;
    void                    read(const GFitsHDU& hdu);
    GFitsHDU*               write(GFits& file,
                                  const std::string& extname = "") const;
    void                    publish(const std::string& name = "") const;
    std::string             print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void              init_members(void);
    void              copy_members(const GSkyMap& map);
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
    double            solidangle(const GSkyDir& dir1, const GSkyDir& dir2,
                                 const GSkyDir& dir3) const;
    double            solidangle(const GSkyDir& dir1, const GSkyDir& dir2,
                                 const GSkyDir& dir3, const GSkyDir& dir4) const;
    bool              overlaps_circle(const GSkyRegionCircle& region) const;
    bool              is_healpix(const GFitsHDU& hdu) const;
    bool              is_wcs(const GFitsHDU& hdu) const;
    GNdarray          smooth_kernel(const std::string& kernel,
                                    const double&      par) const;

    // Private data area
    int               m_num_pixels; //!< Number of pixels (used for pixel allocation)
    int               m_num_maps;   //!< Number of maps (used for pixel allocation)
    int               m_num_x;      //!< Number of pixels in x direction (only 2D)
    int               m_num_y;      //!< Number of pixels in y direction (only 2D)
    std::vector<int>  m_shape;      //!< Shape of the maps
    GSkyProjection*   m_proj;       //!< Pointer to sky projection
    GNdarray          m_pixels;     //!< Skymap pixels

    // Computation cache
    mutable bool      m_hascache;     //!< Cache is valid
    mutable bool      m_contained;    //!< Sky direction is contained in map
    mutable GSkyDir   m_last_dir;     //!< Last sky direction
    mutable GBilinear m_interpol;     //!< Bilinear interpolator
};


/***********************************************************************//**
 * @brief Binary sky map addition
 *
 * @param[in] map Sky map.
 * @return Sky map to which @p map was added.
 *
 * Returns the sum of two sky maps.
 ***************************************************************************/
inline
GSkyMap GSkyMap::operator+(const GSkyMap& map) const
{
    GSkyMap result = *this;
    result        += map;
    return result;
}


/***********************************************************************//**
 * @brief Binary sky map subtraction
 *
 * @param[in] map Sky map.
 * @return Sky map to which @p map was added.
 *
 * Returns the difference of two sky maps.
 ***************************************************************************/
inline
GSkyMap GSkyMap::operator-(const GSkyMap& map) const
{
    GSkyMap result = *this;
    result        -= map;
    return result;
}


/***********************************************************************//**
 * @brief Binary sky map multiplication
 *
 * @param[in] map Sky map.
 * @return Sky map multiplied by @p map.
 *
 * Returns the product of two sky maps.
 ***************************************************************************/
inline
GSkyMap GSkyMap::operator*(const GSkyMap& map) const
{
    GSkyMap result = *this;
    result        *= map;
    return result;
}


/***********************************************************************//**
 * @brief Binary sky map division
 *
 * @param[in] map Sky map.
 * @return Sky map divided by @p map.
 *
 * Returns the ratio of two sky maps.
 ***************************************************************************/
inline
GSkyMap GSkyMap::operator/(const GSkyMap& map) const
{
    GSkyMap result = *this;
    result        /= map;
    return result;
}


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyMap").
 ***************************************************************************/
inline
std::string GSkyMap::classname(void) const
{
    return ("GSkyMap");
}


/***********************************************************************//**
 * @brief Signals if sky map is empty
 *
 * @return True if sky map is empty, false otherwise.
 *
 * Signals if a sky map has no pixels or maps.
 ***************************************************************************/
inline
bool GSkyMap::is_empty(void) const
{
    return ((m_num_pixels == 0) || (m_num_maps == 0));
}


/***********************************************************************//**
 * @brief Returns number of pixels
 *
 * @return Number of pixels in one sky map.
 *
 * Returns the number of pixels in one sky map.
 ***************************************************************************/
inline
const int& GSkyMap::npix(void) const
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
const int& GSkyMap::nx(void) const
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
const int& GSkyMap::ny(void) const
{
    return m_num_y;
}


/***********************************************************************//**
 * @brief Returns number of maps
 *
 * @return Number of maps in the sky map object.
 ***************************************************************************/
inline
const int& GSkyMap::nmaps(void) const
{
    return m_num_maps;
}


/***********************************************************************//**
 * @brief Returns shape of maps
 *
 * @return Shape of maps in the sky map object.
 ***************************************************************************/
inline
const std::vector<int>& GSkyMap::shape(void) const
{
    return m_shape;
}


/***********************************************************************//**
 * @brief Returns dimension of maps
 *
 * @return Number of map dimensions.
 ***************************************************************************/
inline
int GSkyMap::ndim(void) const
{
    return (int)m_shape.size();
}


/***********************************************************************//**
 * @brief Returns pointer to sky projection
 *
 * @return Pointer to sky projection (NULL if no projection is defined).
 ***************************************************************************/
inline
const GSkyProjection* GSkyMap::projection(void) const
{
    return m_proj;
}


/***********************************************************************//**
 * @brief Returns pointer to pixel data
 *
 * @return Pointer to pixel data.
 ***************************************************************************/
inline
const double* GSkyMap::pixels(void) const
{
    return (m_pixels.data());
}

#endif /* GSKYMAP_HPP */
