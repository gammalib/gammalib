/***************************************************************************
 *            GSkymap.hpp  -  Class that implements a sky map              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @brief Sky map class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GSKYMAP_HPP
#define GSKYMAP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GWcs.hpp"
#include "GSkyDir.hpp"
#include "GSkyPixel.hpp"
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
 * @brief GSkymap class interface defintion
 *
 * This class implements a set of skymaps. Skymaps may either be present in
 * the HEALPix format or in any (supported) WCS format. Skymap pixels are
 * stored in a double precision array indexed by (x,y,map), with the x
 * axis being the most rapidely varying axis. 
 * 
 * Skymap pixels may be either accessed via 1D operators (i,map), where i 
 * runs over the (x,y) direction of the map, or via 2D operators (pixel,map),
 * where pixel implements a 2D (x,y) pixel direction. The first operator is
 * the preferred access method for HEALPix maps while the second operator is
 * the preferred access method for WCS maps. Conversion methods between the
 * index or sky pixel and the true physical sky direction are provided by the
 * pix2dir() and dir2pix() methods.
 ***************************************************************************/
class GSkymap : public GBase {

public:
    // Constructors and destructors
    GSkymap(void);
    explicit GSkymap(const std::string& filename);
    explicit GSkymap(const std::string& wcs, const std::string& coords,
                     const int& nside, const std::string& order,
                     const int nmaps = 1);
    explicit GSkymap(const std::string& wcs, const std::string& coords,
                     double const& x, double const& y,
                     double const& dx, double const& dy,
                     const int& nx, const int& ny, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // Operators
    GSkymap& operator= (const GSkymap& map);

    // 1D pixel methods
    double&       operator()(const int& pixel, const int& map = 0);
    const double& operator()(const int& pixel, const int& map = 0) const;
    GSkyDir       pix2dir(const int& pix) const;
    int           dir2pix(const GSkyDir& dir) const;
    double        omega(const int& pix) const;

    // 2D pixel methods
    double&       operator()(const GSkyPixel& pixel, const int& map = 0);
    const double& operator()(const GSkyPixel& pixel, const int& map = 0) const;
    GSkyDir       xy2dir(const GSkyPixel& pix) const;
    GSkyPixel     dir2xy(const GSkyDir& dir) const;
    double        omega(const GSkyPixel& pix) const;

    // Sky direction methods
    double        operator() (const GSkyDir& dir, const int& map = 0) const;

    // Methods
    void          clear(void);
    GSkymap*      clone(void) const;
    void          load(const std::string& filename);
    void          save(const std::string& filename, bool clobber = false) const;
    void          read(const GFitsHDU* hdu);
    void          write(GFits* file) const;
    int           npix(void) const;
    int           nx(void) const;
    int           ny(void) const;
    int           nmaps(void) const;
    int           xy2pix(const GSkyPixel& pix) const;
    GSkyPixel     pix2xy(const int& pix) const;
    GWcs*         wcs(void) const { return m_wcs; }
    void          wcs(const GWcs& wcs);
    double*       pixels(void) const { return m_pixels; }
    bool          isinmap(const GSkyDir& dir) const;
    bool          isinmap(const GSkyPixel& pixel) const;
    std::string   print(void) const;

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
    int     m_num_pixels;   //!< Number of pixels (used for pixel allocation)
    int     m_num_maps;     //!< Number of maps (used for pixel allocation)
    int     m_num_x;        //!< Number of pixels in x direction (only 2D)
    int     m_num_y;        //!< Number of pixels in y direction (only 2D)
    GWcs*   m_wcs;          //!< Pointer to WCS projection
    double* m_pixels;       //!< Pointer to skymap pixels
};

#endif /* GSKYMAP_HPP */
