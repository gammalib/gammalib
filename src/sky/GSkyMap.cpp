/***************************************************************************
 *                       GSkyMap.cpp - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2023 by Juergen Knoedlseder                         *
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
 * @file GSkyMap.cpp
 * @brief Sky map class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GTools.hpp"
#include "GFft.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GSkyDirs.hpp"
#include "GSkyMap.hpp"
#include "GHealpix.hpp"
#include "GWcsRegistry.hpp"
#include "GWcs.hpp"
#include "GMatrix.hpp"
#include "GVector.hpp"
#include "GVOClient.hpp"
#include "GSkyRegion.hpp"
#include "GSkyRegionCircle.hpp"
#include "GSkyRegions.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT_HPX                "GSkyMap::GSkyMap(std::string&, int&,"\
                                                       " std::string&, int&)"
#define G_CONSTRUCT_MAP        "GSkyMap::GSkyMap(std::string&, std::string&,"\
                      " double&, double&, double& double&, int&, int&, int&)"
#define G_NMAPS                                        "GSkyMap::nmaps(int&)"
#define G_SHAPE                           "GSkyMap::shape(std::vector<int>&)"
#define G_OP_UNARY_ADD                        "GSkyMap::operator+=(GSkyMap&)"
#define G_OP_UNARY_SUB                        "GSkyMap::operator-=(GSkyMap&)"
#define G_OP_UNARY_MUL                        "GSkyMap::operator-=(GSkyMap&)"
#define G_OP_UNARY_DIV                        "GSkyMap::operator/=(GSkyMap&)"
#define G_OP_UNARY_DIV2                        "GSkyMap::operator/=(double&)"
#define G_OP_ACCESS_1D                        "GSkyMap::operator(int&, int&)"
#define G_OP_ACCESS_2D                  "GSkyMap::operator(GSkyPixel&, int&)"
#define G_OP_VALUE                        "GSkyMap::operator(GSkyDir&, int&)"
#define G_INX2DIR                                    "GSkyMap::inx2dir(int&)"
#define G_PIX2DIR                              "GSkyMap::pix2dir(GSkyPixel&)"
#define G_DIR2INX                                "GSkyMap::dir2inx(GSkyDir&)"
#define G_DIR2PIX                                "GSkyMap::dir2pix(GSkyDir&)"
#define G_FLUX1                                         "GSkyMap::flux(int&)"
#define G_FLUX2                                   "GSkyMap::flux(GSkyPixel&)"
#define G_FLUX3                            "GSkyMap::flux(GSkyRegion&, int&)"
#define G_FLUX4                           "GSkyMap::flux(GSkyRegions&, int&)"
#define G_SOLIDANGLE1                             "GSkyMap::solidangle(int&)"
#define G_SOLIDANGLE2                       "GSkyMap::solidangle(GSkyPixel&)"
#define G_OVERLAPS                           "GSkyMap::overlaps(GSkyRegion&)"
#define G_EXTRACT                              "GSkyMap::extract(int&, int&)"
#define G_EXTRACT_INT              "GSkyMap::extract(int&, int&, int&, int&)"
#define G_EXTRACT_REG                        "GSkyMap::extract(GSkyRegions&)"
#define G_READ                                     "GSkyMap::read(GFitsHDU&)"
#define G_SET_WCS     "GSkyMap::set_wcs(std::string&, std::string&, double&,"\
                              " double&, double&, double&, double&, double&)"
#define G_READ_HEALPIX                   "GSkyMap::read_healpix(GFitsTable*)"
#define G_READ_WCS                           "GSkyMap::read_wcs(GFitsImage*)"
#define G_ALLOC_WCS                         "GSkyMap::alloc_wcs(GFitsImage*)"
#define G_SQRT                                       "GSkyMap sqrt(GSkyMap&)"
#define G_LOG                                         "GSkyMap log(GSkyMap&)"
#define G_LOG10                                     "GSkyMap log10(GSkyMap&)"
#define G_CONVOLUTION_KERNEL     "GSkyMap::convolution_kernel(std::string&, "\
                                                            "double&, bool&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_READ_HEALPIX_DEBUG                          // Debug read_healpix
//#define G_READ_WCS_DEBUG                                  // Debug read_wcs

/* __ Prototype __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                          Constructors/destructors                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty sky map.
 ***************************************************************************/
GSkyMap::GSkyMap(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS file constructor
 *
 * @param[in] filename FITS file name.
 *
 * Constructs a sky map by loading data from a FITS file. See the load()
 * method for more information about the FITS formats that are supported.
 ***************************************************************************/
GSkyMap::GSkyMap(const GFilename& filename)
{
    // Initialise class members for clean destruction
    init_members();

    // Load skymap
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief FITS HDU constructor
 *
 * @param[in] hdu FITS HDU.
 *
 * Constructs a sky map by fetching data from a FITS HDU. See the read()
 * method for more information.
 ***************************************************************************/
GSkyMap::GSkyMap(const GFitsHDU& hdu)
{
    // Initialise class members for clean destruction
    init_members();

    // Read skymap from HDU
    read(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix sky map constructor
 *
 * @param[in] coords Coordinate System (CEL or GAL).
 * @param[in] nside Nside parameter.
 * @param[in] order Pixel ordering (RING or NEST).
 * @param[in] nmaps Number of maps in set.
 *
 * @exception GException::invalid_argument
 *            Invalid @p nmaps parameter.
 *
 * Constructs @p nmaps identical all sky maps in Healpix pixelisation. All
 * pixels of the sky maps will be initialised to values of zero.
 ***************************************************************************/
GSkyMap::GSkyMap(const std::string& coords,
                 const int&         nside,
                 const std::string& order,
                 const int&         nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    // Check if nmaps parameter is >0
    if (nmaps < 1) {
        std::string msg = "Number of maps "+gammalib::str(nmaps)+" is smaller "
                          "than 1. Please specify at least one map.";
        throw GException::invalid_argument(G_CONSTRUCT_HPX, msg);
    }

    // Allocate Healpix projection
    GHealpix* projection = new GHealpix(nside, order, coords);
    m_proj               = projection;

    // Set number of pixels, number of maps and shape of maps
    m_num_pixels = projection->npix();
    m_num_maps   = nmaps;
    m_shape.push_back(nmaps);

    // Allocate pixels
    m_pixels = GNdarray(m_num_pixels, m_num_maps);

    // Return
    return;
}


/***********************************************************************//**
 * @brief WCS sky map constructor
 *
 * @param[in] wcs World Coordinate System.
 * @param[in] coords Coordinate System (CEL or GAL).
 * @param[in] x X coordinate of sky map centre (deg).
 * @param[in] y Y coordinate of sky map centre (deg).
 * @param[in] dx Pixel size in x direction at centre (deg/pixel).
 * @param[in] dy Pixel size in y direction at centre (deg/pixel).
 * @param[in] nx Number of pixels in x direction.
 * @param[in] ny Number of pixels in y direction.
 * @param[in] nmaps Number of maps in set (default=1).
 *
 * @exception GException::invalid_argument
 *            Invalid number of pixels or maps.
 *
 * Constructs @p nmaps identical all sky maps in World Coordinate System
 * projection. All pixels of the sky maps will be initialised to values of
 * zero.
 ***************************************************************************/
GSkyMap::GSkyMap(const std::string& wcs,
                 const std::string& coords,
                 const double&      x,
                 const double&      y,
                 const double&      dx,
                 const double&      dy,
                 const int&         nx,
                 const int&         ny,
                 const int&         nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    // Check parameters
    if (nx < 1) {
        std::string msg = "Number of pixels "+gammalib::str(nx)+" in x "
                          "direction is smaller than 1. Please specify at "
                          "least one pixel in x direction.";
        throw GException::invalid_argument(G_CONSTRUCT_MAP, msg);
    }
    if (ny < 1) {
        std::string msg = "Number of pixels "+gammalib::str(ny)+" in y "
                          "direction is smaller than 1. Please specify at "
                          "least one pixel in y direction.";
        throw GException::invalid_argument(G_CONSTRUCT_MAP, msg);
    }
    if (nmaps < 1) {
        std::string msg = "Number of maps "+gammalib::str(nmaps)+" is smaller "
                          "than 1. Please specify at least one map.";
        throw GException::invalid_argument(G_CONSTRUCT_MAP, msg);
    }

    // Set WCS
    double  crval1 = x;
    double  crval2 = y;
    double  crpix1 = double(nx+1)/2.0;
    double  crpix2 = double(ny+1)/2.0;
    double  cdelt1 = dx;
    double  cdelt2 = dy;
    set_wcs(wcs, coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);

    // Set number of pixels, number of maps and shape of maps
    m_num_x      = nx;
    m_num_y      = ny;
    m_num_pixels = m_num_x * m_num_y;
    m_num_maps   = nmaps;
    m_shape.push_back(nmaps);

    // Allocate pixels
    m_pixels = GNdarray(nx, ny, nmaps);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] map Sky map.
 *
 * Constructs sky maps by copying data from another sky map object.
 ***************************************************************************/
GSkyMap::GSkyMap(const GSkyMap& map)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(map);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyMap::~GSkyMap(void)
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
 * @param[in] map Sky map.
 * @return Sky map.
 *
 * Assigns one sky map to another.
 ***************************************************************************/
GSkyMap& GSkyMap::operator=(const GSkyMap& map)
{
    // Execute only if object is not identical
    if (this != &map) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(map);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Value setting operator
 *
 * @param[in] value Value.
 * @return Sky map.
 *
 * Sets all sky map pixels to the specified @p value.
 ***************************************************************************/
GSkyMap& GSkyMap::operator=(const double& value)
{
    // Get number of pixels
    int num = m_pixels.size();

    // Loop over all pixels
    for (int i = 0; i < num; ++i) {
        m_pixels(i) = value;
    }

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Map addition operator
 *
 * @param[in] map Sky map.
 * @return Sky map.
 *
 * @exception GException::invalid_value
 *            Mismatch between number of maps in skymap object.
 *
 * Adds the content of @p map to the skymap. The operator only works on sky
 * maps with an identical number of maps.
 *
 * If the sky maps have the same definitions, the operator will add the
 * values of the sky map pixels. Otherwise, the map values are added by
 * bi-linearily interpolating the values in the source sky map, allowing thus
 * for a reprojection of sky map values.
 ***************************************************************************/
GSkyMap& GSkyMap::operator+=(const GSkyMap& map)
{
    // Check if number of layers are identical
    if (map.nmaps() != nmaps()) {
        std::string msg = "Mismatch of number of maps in skymap object"
                          " ("+gammalib::str(nmaps())+" maps in destination"
                          " map, "+gammalib::str(map.nmaps())+" in source"
                          " map).";
        throw GException::invalid_value(G_OP_UNARY_ADD, msg);
    }

    // If maps are identical then add the pixels
    if (is_same(map)) {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Add pixel
                (*this)(index, layer) += map(index, layer);

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    }

    // ... otherwise use interpolation method and add map
    else {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Get sky direction of actual pixel
            GSkyDir dir = inx2dir(index);

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Add value
                (*this)(index, layer) += map(dir, layer);

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    } // endelse: interpolation method

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Value addition operator
 *
 * @param[in] value Value.
 * @return Sky map.
 *
 * Add @p value to all sky map pixels.
 ***************************************************************************/
GSkyMap& GSkyMap::operator+=(const double& value)
{
    // Add value to all pixels
    m_pixels += value;

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Map subtraction operator
 *
 * @param[in] map Sky map.
 * @return Sky map.
 *
 * @exception GException::invalid_value
 *            Mismatch between number of maps in skymap object.
 *
 * Subtracts the content of @p map from the skymap. The operator only works
 * on sky maps with an identical number of layers.
 *
 * If the sky maps have the same definitions, the operator will subtract the
 * values of the sky map pixels. Otherwise, the map values are subtracted by
 * bi-linearily interpolating the values in the source sky map, allowing thus
 * for a reprojection of sky map values.
 ***************************************************************************/
GSkyMap& GSkyMap::operator-=(const GSkyMap& map)
{
    // Check if number of layers are identical
    if (map.nmaps() != nmaps()) {
        std::string msg = "Mismatch of number of maps in skymap object"
                          " ("+gammalib::str(nmaps())+" maps in destination"
                          " map, "+gammalib::str(map.nmaps())+" in source"
                          " map).";
        throw GException::invalid_value(G_OP_UNARY_SUB, msg);
    }

    // If maps are identical then subtract the pixels
    if (is_same(map)) {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Subtract pixel
                (*this)(index, layer) -= map(index, layer);

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    }

    // ... otherwise use interpolation method and subtract map
    else {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Get sky direction of actual pixel
            GSkyDir dir = inx2dir(index);

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Subtract value
                (*this)(index, layer) -= map(dir, layer);

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    } // endelse: interpolation method

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Value subtraction operator
 *
 * @param[in] value Value.
 * @return Sky map.
 *
 * Subtracts @p value from all sky map pixels.
 ***************************************************************************/
GSkyMap& GSkyMap::operator-=(const double& value)
{
    // Subtract value from all pixels
    m_pixels -= value;

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Multiplication operator
 *
 * @param[in] map Sky map.
 * @return Sky map.
 *
 * @exception GException::invalid_value
 *            Mismatch between number of maps in skymap object.
 *
 * Multiplies the content of @p map from the skymap. The operator only works
 * on sky maps with an identical number of layers.
 *
 * If the sky maps have the same definitions, the operator will multiply the
 * values of the sky map pixels. Otherwise, the map values are multiplied by
 * bi-linearily interpolating the values in the source sky map, allowing thus
 * for a reprojection of sky map values.
 ***************************************************************************/
GSkyMap& GSkyMap::operator*=(const GSkyMap& map)
{
    // Check if number of layers are identical
    if (map.nmaps() != nmaps()) {
        std::string msg = "Mismatch of number of maps in skymap object"
                          " ("+gammalib::str(nmaps())+" maps in destination"
                          " map, "+gammalib::str(map.nmaps())+" in source"
                          " map).";
        throw GException::invalid_value(G_OP_UNARY_MUL, msg);
    }

    // If maps are identical then multiply the pixels
    if (is_same(map)) {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Multiply pixel
                (*this)(index, layer) *= map(index, layer);

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    }

    // ... otherwise use interpolation method and multiply map
    else {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Get sky direction of actual pixel
            GSkyDir dir = inx2dir(index);

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Multiply value
                (*this)(index, layer) *= map(dir, layer);

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    } // endelse: interpolation method

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Sky map scaling operator
 *
 * @param[in] factor Scale factor.
 * @return Sky map.
 *
 * Multiplies all pixels of the sky map by the given scale @p factor.
 ***************************************************************************/
GSkyMap& GSkyMap::operator*=(const double& factor)
{
    // Scale all pixels
    m_pixels *= factor;

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Division operator
 *
 * @param[in] map Sky map.
 * @return Sky map.
 *
 * @exception GException::invalid_value
 *            Mismatch between number of maps in skymap object.
 *
 * Divides the content of the actual skymap by the skymap @p map. The operator
 * only works on sky maps with an identical number of layers.
 *
 * If the sky maps have the same definitions, the operator will divide the
 * values of the sky map pixels. Otherwise, the map values are divided by
 * bi-linearily interpolating the values in the source sky map, allowing thus
 * for a reprojection of sky map values.
 *
 * On return, all pixels in @p map that are zero are silently set to zero in
 * the skymap, hence no division by zero exception will be raised in case
 * that zero values are encountered.
 ***************************************************************************/
GSkyMap& GSkyMap::operator/=(const GSkyMap& map)
{
    // Check if number of layers are identical
    if (map.nmaps() != nmaps()) {
        std::string msg = "Mismatch of number of maps in skymap object"
                          " ("+gammalib::str(nmaps())+" maps in destination"
                          " map, "+gammalib::str(map.nmaps())+" in source"
                          " map).";
        throw GException::invalid_value(G_OP_UNARY_DIV, msg);
    }

    // If maps are identical then divide the pixels
    if (is_same(map)) {

        // Loop over all pixels of sky map
        for (int index = 0; index < npix(); ++index) {

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Get pixel value of the map by which the division will be
                // done.
                double value = map(index, layer);

                // Divide destination pixel by value. In case of a division
                // by zero set the destination pixel also to zero.
                if (value == 0.0) {
                    (*this)(index, layer) = 0.0;
                }
                else {
                    (*this)(index, layer) /= value;
                }

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    }

    // ... otherwise use interpolation method and divide map
    else {

        // Loop over all pixels of destination sky map
        for (int index = 0; index < npix(); ++index) {

            // Get sky direction of actual pixel
            GSkyDir dir = inx2dir(index);

            // Loop over all layers
            for (int layer = 0; layer < nmaps(); ++layer) {

                // Get map value of the map by which the division will be
                // done. The map is accessed using the sky direction, hence
                // interpolating properly the map. If we're outside the map,
                // zero is returned.
                double value = map(dir, layer);

                // Divide destination pixel by value. In case of a division
                // by zero set the destination pixel also to zero.
                if (value == 0.0) {
                    (*this)(index, layer) = 0.0;
                }
                else {
                    (*this)(index, layer) /= value;
                }

            } // endfor: looped over all layers

        } // endfor: looped over all pixels

    } // endelse: interpolation method

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Sky map division operator
 *
 * @param[in] factor Scale factor.
 * @return Sky map.
 *
 * @exception GException::invalid_argument
 *            Division by zero error.
 *
 * Divides all pixels of the sky map by the given @p factor.
 ***************************************************************************/
GSkyMap& GSkyMap::operator/=(const double& factor)
{
    // Check for division by zero
    if (factor == 0.0) {
        std::string msg = "Trying to divide sky map pixels by zero.";
        throw GException::invalid_argument(G_OP_UNARY_DIV2, msg);
    }

    // Divide all pixels by factor
    m_pixels /= factor;

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Pixel index access operator
 *
 * @param[in] index Pixel index [0,...,npix()-1].
 * @param[in] map Map index [0,...,nmaps()-1].
 * @return Sky map pixel value.
 *
 * @exception GException::out_of_range
 *            Pixel index and/or map index are outside valid range.
 *
 * Access sky map pixel by its index, where the most quickly varying axis is
 * the x axis of the map.
 ***************************************************************************/
double& GSkyMap::operator()(const int& index, const int& map)
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num_pixels) {
        throw GException::out_of_range(G_OP_ACCESS_1D,
                                       "Sky map pixel index",
                                       index, m_num_pixels);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_1D,
                                       "Sky map map index",
                                       map, m_num_maps);
    }
    #endif

    // Return reference to pixel value
    return m_pixels(index+m_num_pixels*map);
}


/***********************************************************************//**
 * @brief Pixel index access operator (const variant)
 *
 * @param[in] index Pixel index [0,...,npix()-1].
 * @param[in] map Map index [0,...,nmaps()-1].
 *
 * @exception GException::out_of_range
 *            Pixel index and/or map index are outside valid range.
 *
 * Access sky map pixel by its index, where the most quickly varying axis is
 * the x axis of the map.
 ***************************************************************************/
const double& GSkyMap::operator()(const int& index, const int& map) const
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_num_pixels) {
        throw GException::out_of_range(G_OP_ACCESS_1D,
                                       "Sky map pixel index",
                                       index, m_num_pixels);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_1D,
                                       "Sky map map index",
                                       map, m_num_maps);
    }
    #endif

    // Return reference to pixel value
    return m_pixels(index+m_num_pixels*map);
}


/***********************************************************************//**
 * @brief Sky map pixel access operator
 *
 * @param[in] pixel Sky map pixel.
 * @param[in] map Sky map index [0,...,nmaps()[.
 *
 * @exception GException::invalid_argument
 *            Sky pixel is not contained in sky map.
 * @exception GException::out_of_range
 *            Sky map index is out of range.
 *
 * Access sky map pixel by its 2D index (x,y) that is implemented by the
 * GSkyPixel class.
 ***************************************************************************/
double& GSkyMap::operator()(const GSkyPixel& pixel, const int& map)
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (!contains(pixel)) {
        std::string msg = "Sky pixel ("+gammalib::str(pixel.x())+","+
                          gammalib::str(pixel.y())+") is not contained in "
                          "sky map comprised of ("+gammalib::str(m_num_x)+","+
                          gammalib::str(m_num_y)+") pixels. Please specify "
                          "a valid sky pixel.";
        throw GException::invalid_argument(G_OP_ACCESS_2D, msg);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_2D, "Sky map index",
                                       map, m_num_maps);
    }
    #endif

    // Get pixel index
    int index = pix2inx(pixel);

    // Return reference to pixel value
    return m_pixels(index+m_num_pixels*map);
}


/***********************************************************************//**
 * @brief Sky map pixel access operator
 *
 * @param[in] pixel Sky map pixel.
 * @param[in] map Sky map index [0,...,nmaps()[.
 *
 * @exception GException::invalid_argument
 *            Sky pixel is not contained in sky map.
 * @exception GException::out_of_range
 *            Sky map index is out of range.
 *
 * Access sky map pixel by its 2D index (x,y) that is implemented by the
 * GSkyPixel class.
 ***************************************************************************/
const double& GSkyMap::operator()(const GSkyPixel& pixel, const int& map) const
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (!contains(pixel)) {
        std::string msg = "Sky pixel ("+gammalib::str(pixel.x())+","+
                          gammalib::str(pixel.y())+") is not contained in "
                          "sky map comprised of ("+gammalib::str(m_num_x)+","+
                          gammalib::str(m_num_y)+") pixels. Please specify "
                          "a valid sky pixel.";
        throw GException::invalid_argument(G_OP_ACCESS_2D, msg);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_2D, "Sky map index",
                                       map, m_num_maps);
    }
    #endif

    // Get pixel index
    int index = pix2inx(pixel);

    // Return reference to pixel value
    return m_pixels(index+m_num_pixels*map);
}


/***********************************************************************//**
 * @brief Return interpolated skymap value for sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] map Map index [0,...,nmaps()-1].
 * @return Sky intensity.
 *
 * @exception GException::out_of_range
 *            Map index lies outside valid range.
 *
 * Returns the skymap value for a given sky direction, obtained by bi-linear
 * interpolation of the neighbouring pixels. If the sky direction falls
 * outside the area covered by the skymap, a value of 0 is returned.
 *
 * The method implements a computation cache that avoids recomputation of
 * interpolation indices and weights in case that the same sky direction
 * is requested several times. This speeds up use cases where skymap values
 * for various map indices have to be returned for the same sky direction.
 ***************************************************************************/
double GSkyMap::operator()(const GSkyDir& dir, const int& map) const
{
    // Throw an error if the map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_VALUE,
                                       "Sky map map index",
                                       map, m_num_maps);
    }
    #endif

    // Initialise intensity
    double intensity = 0.0;

    // If the pre-computation cache has not been set or if the sky direction
    // has changed then compute the bilinear interpolator
    if (!m_hascache || (dir != m_last_dir)) {

        // Perform computation for HealPix map
        if (m_proj->size() == 1)  {

            // Get Healpix interpolator
            m_interpol = static_cast<GHealpix*>(m_proj)->interpolator(dir);

            // Set containment flag and store actual sky direction as last
            // direction
            m_hascache  = true;
            m_contained = true;
            m_last_dir  = dir;

        } // endif: we had a Healpix projection

        // ... otherwise perform computation for WCS map
        else {

            // Determine sky pixel. At this point an exception may occur
            // in case that the pixel cannot be represented by the
            // relevant projection. We catch this exception here
            try {

                // Determine sky pixel
                GSkyPixel pixel = dir2pix(dir);

                // Set containment flag and store actual sky direction
                // as last direction
                m_hascache  = true;
                m_contained = contains(pixel);
                m_last_dir  = dir;

                // Continue only if pixel is within the map
                if (m_contained) {

                    // Signal that we have a wrap around in the x axis
                    double x_size = std::abs(m_num_x * static_cast<GWcs*>(m_proj)->cdelt(0));
                    bool   x_wrap = (x_size > 359.99);

                    // Get pixel index that is to the left-top of the actual
                    // pixel. We take care of the special case of negative pixel
                    // indices which arise if we are in the first column or
                    // first row of the map.
                    int inx_x = int(pixel.x());
                    int inx_y = int(pixel.y());
                    if (pixel.x() < 0.0) {
                        inx_x--;
                    }
                    if (pixel.y() < 0.0) {
                        inx_y--;
                    }

                    // Set left and right indices for interpolation. The left
                    // index needs to be non-negative and the right index needs
                    // to be not larger that the number of pixels. We treat
                    // here also the special case of wrap around in the x
                    // axis that may occur if we have an allsky map.
                    int inx_x_left  = inx_x;
                    int inx_y_left  = inx_y;
                    int inx_x_right = inx_x_left + 1;
                    int inx_y_right = inx_y_left + 1;
                    if (inx_x_left < 0) {
                        if (x_wrap) {
                            inx_x_left += m_num_x;
                        }
                        else {
                            inx_x_left  = 0;
                        }
                    }
                    if (inx_x_right >= m_num_x) {
                        if (x_wrap) {
                            inx_x_right -= m_num_x;
                        }
                        else {
                            inx_x_right = m_num_x - 1;
                        }
                    }
                    if (inx_y_left < 0) {
                        inx_y_left  = 0;
                    }
                    if (inx_y_right >= m_num_y) {
                        inx_y_right = m_num_y - 1;
                    }

                    // Set weighting factors for interpolation
                    double wgt_x_right = (pixel.x() - inx_x);
                    double wgt_x_left  = 1.0 - wgt_x_right;
                    double wgt_y_right = (pixel.y() - inx_y);
                    double wgt_y_left  = 1.0 - wgt_y_right;

                    // Compute skymap pixel indices for bi-linear interpolation
                    m_interpol.index1() = inx_x_left  + inx_y_left  * m_num_x;
                    m_interpol.index2() = inx_x_left  + inx_y_right * m_num_x;
                    m_interpol.index3() = inx_x_right + inx_y_left  * m_num_x;
                    m_interpol.index4() = inx_x_right + inx_y_right * m_num_x;

                    // Compute weighting factors for bi-linear interpolation
                    m_interpol.weight1() = wgt_x_left  * wgt_y_left;
                    m_interpol.weight2() = wgt_x_left  * wgt_y_right;
                    m_interpol.weight3() = wgt_x_right * wgt_y_left;
                    m_interpol.weight4() = wgt_x_right * wgt_y_right;

                } // endif: pixel was contained in map

            } // endtry: pixel computation was successful
            catch (GException::invalid_argument) {
                m_contained = false;
            }

        } // endelse: we had a WCS map

    } // endif: computation of the bi-linear interpolater was required

    // Compute the interpolated intensity if the pixel is contained in
    // the map
    if (m_contained) {

        // Compute map offset
        int offset = m_num_pixels * map;

        // Compute interpolated skymap value
        intensity = m_interpol(m_pixels.data()+offset);

    } // endif: direction was contained in map

    // Return intensity
    return intensity;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

 /***********************************************************************//**
 * @brief Clear instance.
 *
 * Resets the sky map to the initial state.
 ***************************************************************************/
void GSkyMap::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone sky map
 *
 * @return Pointer to deep copy of sky map.
 ***************************************************************************/
GSkyMap* GSkyMap::clone(void) const
{
    return new GSkyMap(*this);
}


/***********************************************************************//**
 * @brief Set number of maps
 *
 * @param[in] nmaps Number of maps.
 *
 * @exception GException::invalid_argument
 *            Invalid number of maps specified.
 *
 * Redefines the number of maps in an GSkyMap object. If the number of maps
 * is increased with respect to the existing number, additional maps with
 * pixel values of zero are append to the object. Existing map pixel values
 * are kept. If the number of maps is decreased with respect to the existing
 * number, the excedent maps are dropped. The remaining map pixel values are
 * kept.
 ***************************************************************************/
void GSkyMap::nmaps(const int& nmaps)
{
    // Throw an exception if less than 1 map is required
    if (nmaps < 1) {
        std::string msg = "At least one map is required in an GSkyMap object."
                          " Please specify an argument >=1.";
        throw GException::invalid_argument(G_NMAPS, msg);
    }

    // If the map has pixels and if the number of maps changes then set a new
    // number of maps. Copy over any existing information.
    if (m_num_pixels > 0 && nmaps != m_num_maps) {

        // Allocate new array
        GNdarray new_pixels(m_num_pixels, nmaps);

        // Copy over existing pixels
        int num_copy = (nmaps > m_num_maps) ? m_num_maps : nmaps;
        num_copy    *= m_num_pixels;
        for (int i = 0; i < num_copy; ++i) {
            new_pixels(i) = m_pixels(i);
        }

        // Set pixels to new array
        m_pixels = new_pixels;

        // Set number of maps
        m_num_maps = nmaps;

        // Set shape of maps to linear array since we do not know how to
        // arrange the new maps
        m_shape.clear();
        m_shape.push_back(nmaps);

    } // endif: map had pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set one-dimensional shape of maps
 *
 * @param[in] s1 Axis length of first dimension.
 ***************************************************************************/
void GSkyMap::shape(const int& s1)
{
    // Initialise vector
    std::vector<int> shape;

    // Set vector elements
    shape.push_back(s1);

    // Call vector method
    this->shape(shape);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set two-dimensional shape of maps
 *
 * @param[in] s1 Axis length of first dimension.
 * @param[in] s2 Axis length of second dimension.
 ***************************************************************************/
void GSkyMap::shape(const int& s1, const int& s2)
{
    // Initialise vector
    std::vector<int> shape;

    // Set vector elements
    shape.push_back(s1);
    shape.push_back(s2);

    // Call vector method
    this->shape(shape);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set three-dimensional shape of maps
 *
 * @param[in] s1 Axis length of first dimension.
 * @param[in] s2 Axis length of second dimension.
 * @param[in] s3 Axis length of second dimension.
 ***************************************************************************/
void GSkyMap::shape(const int& s1, const int& s2, const int& s3)
{
    // Initialise vector
    std::vector<int> shape;

    // Set vector elements
    shape.push_back(s1);
    shape.push_back(s2);
    shape.push_back(s3);

    // Call vector method
    this->shape(shape);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set shape of maps
 *
 * @param[in] shape Shape vector.
 *
 * @exception GException::invalid_argument
 *            Invalid shape factorisation specified.
 *
 * Defines a shape for the maps in the object. The shape specifies how the
 * maps are arranged in a n-dimensional array.
 ***************************************************************************/
void GSkyMap::shape(const std::vector<int>& shape)
{
    // Computes the resulting number of maps
    int nmaps = 0;
    if (shape.size() > 0) {
        nmaps = 1;
        for (int i = 0; i < shape.size(); ++i) {
            nmaps *= shape[i];
        }
    }

    // Throw an exception if resulting number of maps is not equal to the
    // existing number of maps
    if (nmaps != m_num_maps) {
        std::string msg = "Sky map shaping requires "+gammalib::str(nmaps)+
                          " maps while "+gammalib::str(m_num_maps)+" maps "
                          "are available. Please specify an appropriate "
                          "factorisation or modify the number of available "
                          "maps.";
        throw GException::invalid_argument(G_SHAPE, msg);
    }

    // Set shape
    m_shape = shape;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Converts pixel index into sky map pixel
 *
 * @param[in] index Pixel index [0,...,npix()-1].
 * @return Sky map pixel.
 *
 * Converts the pixel @p index into a sky map pixel (GSkyPixel).The dimension
 * of GSkyPixel will be identical to the dimension of the sky map (i.e. a 1D
 * sky map leads to a 1D GSkyPixel object, a 2D sky map leads to a 2D
 * GSkyPixel object.
 ***************************************************************************/
GSkyPixel GSkyMap::inx2pix(const int& index) const
{
    // Initialise sky map pixel
    GSkyPixel pixel;

    // Get x and y indices
    if (m_num_x != 0) { //!< 2D sky map
        pixel.x(double(index % m_num_x));
        pixel.y(double(index / m_num_x));
    }
    else {              //!< 1D sky map
        pixel.index(index);
    }

    // Return pixel
    return pixel;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] index Pixel index [0,...,npix()-1].
 * @return Sky direction.
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns sky direction for a given pixel index.
 ***************************************************************************/
GSkyDir GSkyMap::inx2dir(const int& index) const
{
    // Throw error if sky projection is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_INX2DIR, msg);
    }

    // Determine sky direction from pixel index.
    GSkyDir dir = m_proj->pix2dir(inx2pix(index));

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pixel Sky map pixel.
 * @return Sky direction.
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 * @exception GException::invalid_argument
 *            2D sky map pixel used to access 1D projection.
 *
 * Returns sky direction for a given sky map @p pixel.
 ***************************************************************************/
GSkyDir GSkyMap::pix2dir(const GSkyPixel& pixel) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_PIX2DIR, msg);
    }

    // Initialise sky direction
    GSkyDir dir;

    // If pixel size matches the projection size then perform a straight
    // forward conversion
    if (m_proj->size() == pixel.size()) {
        dir = m_proj->pix2dir(pixel);
    }

    // ... otherwise, if we have a 2D projection but a 1D pixel then
    // interpret the pixel as the linear index in the pixel array
    else if (m_proj->size() == 2) {
        dir = m_proj->pix2dir(GSkyPixel(inx2pix(int(pixel))));
    }

    // ... otherwise we have a 1D projection but a 2D pixel. There is
    // no unambiguous way to handle this case, hence we throw an exception
    else {
        std::string msg = "A 2-dimensional sky map pixel "+pixel.print()+
                          " is used to determine the sky direction for"
                          " the 1-dimensional sky projection \""+
                          m_proj->name()+"\"\n"
                          "Please specify a 1-dimensional sky map pixel.";
        throw GException::invalid_argument(G_PIX2DIR, msg);
    }

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Converts sky map pixel into pixel index
 *
 * @param[in] pixel Sky map pixel.
 * @return Pixel index [0,...,npix()-1].
 *
 * Converts a sky map @p pixel into the pixel index.
 ***************************************************************************/
int GSkyMap::pix2inx(const GSkyPixel& pixel) const
{
    // Initialise pixel index
    int index = 0;

    // Handle 1D sky map pixel
    if (pixel.is_1D()) {
        index = int(pixel);
    }

    // Handle 2D sky map pixel
    else if (pixel.is_2D()) {

        // Get x and y indices by rounding the (x,y) values
        int ix = int(pixel.x()+0.5);
        int iy = int(pixel.y()+0.5);

        // Set index
        index = ix + iy * m_num_x;

    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Returns pixel index for a given sky direction
 *
 * @param[in] dir Sky direction.
 * @return Pixel index [0,...,npix()-1].
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns sky map pixel index for a given sky direction.
 ***************************************************************************/
int GSkyMap::dir2inx(const GSkyDir& dir) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_DIR2INX, msg);
    }

    // Determine pixel index for a given sky direction
    int index = pix2inx(m_proj->dir2pix(dir));

    // Return pixel index
    return index;
}


/***********************************************************************//**
 * @brief Returns sky map pixel for a given sky direction
 *
 * @param[in] dir Sky direction.
 * @return Sky map pixel.
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns sky map pixel for a given sky direction.
 ***************************************************************************/
GSkyPixel GSkyMap::dir2pix(const GSkyDir& dir) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_DIR2PIX, msg);
    }

    // Determine pixel for a given sky direction
    GSkyPixel pixel = m_proj->dir2pix(dir);

    // Return pixel
    return pixel;
}


/***********************************************************************//**
 * @brief Returns array with total number of counts for count maps
 *
 * @return Array of pixel counts sums.
 *
 * For a set of n count maps, the method returns a 1D array with n entries
 * that correspond to the sum of the pixel counts in each map.
 ***************************************************************************/
GNdarray GSkyMap::counts(void) const
{
    // Initialise output array
    GNdarray sum_array = GNdarray(m_num_maps);

    // Loop over maps and store sums in output array
    for (int i = 0; i < m_num_maps; ++i) {

        // Initialise sum for this map
        double sum = 0.0;

        // Compute map sum
        for (int k = 0; k < m_num_pixels; ++k) {
            sum += (*this)(k,i);
        }

        // Store map sum
        sum_array(i) = sum;

    } // endfor: looped over all maps

    // Return sum array
    return sum_array;
}


/***********************************************************************//**
 * @brief Returns flux in pixel
 *
 * @param[in] index Pixel index [0,...,npix()-1].
 * @param[in] map Map index [0,...,nmaps()-1].
 * @return Flux in pixel.
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns the flux in the pixel with the specified @p index. The flux is
 * computed by integrating the intensity over the solid angle subtended by
 * the pixel.
 * Integration is done by dividing the pixel into wedges, computing the
 * intensitites at the cornes of each wedge, averaging these intensities,
 * and multiplying it with the solid angle of each wedge. This provides an
 * approximation of the true pixel flux which is accurate to better than
 * about 5%.
 * 
 * For a HealPix pixelisation, the pixel is divided into 12 wedges. For a
 * WCS pixelisation, the pixel is divided into 8 wedges.
 *
 * @warning
 * This method only returns correct values once the skymap is completely
 * setup with values. Do not use this method during a setup operation as
 * the method uses neighboring pixels for interpolation. If the map is not
 * completely setup the neighboring pixels may still be empty, hence the
 * flux interpolation will be wrong.
 ***************************************************************************/
double GSkyMap::flux(const int& index, const int& map) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_FLUX1, msg);
    }

    // Determine flux from pixel index.
    double flux = this->flux(inx2pix(index), map);

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns flux in pixel
 *
 * @param[in] pixel Sky map pixel.
 * @param[in] map Map index [0,...,nmaps()-1].
 * @return Flux in pixel (ph/cm2/s).
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns the flux in the specified sky map @p pixel. The flux is computed
 * by integrating the intensity over the solid angle subtended by the pixel.
 * Integration is done by dividing the pixel into wedges, computing the
 * intensitites at the cornes of each wedge, averaging these intensities,
 * and multiplying it with the solid angle of each wedge. This provides an
 * approximation of the true pixel flux which is accurate to better than
 * about 5%.
 * 
 * For a HealPix pixelisation, the pixel is divided into 12 wedges. For a
 * WCS pixelisation, the pixel is divided into 8 wedges.
 *
 * @warning
 * This method only returns correct values once the skymap is completely
 * setup with values. Do not use this method during a setup operation as
 * the method uses neighboring pixels for interpolation. If the map is not
 * completely setup the neighboring pixels may still be empty, hence the
 * flux interpolation will be wrong.
 ***************************************************************************/
double GSkyMap::flux(const GSkyPixel& pixel, const int& map) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_FLUX2, msg);
    }

    // Initialise flux
    double flux = 0.0;

    // Perform flux computation for HealPix map. The pixel is divided into
    // 12 wedges and the flux is integrated in each wedge by computing the
    // average of the corner intensities multiplied by the solid angle of
    // the wedge. This leads to a precision of the order of 1-2% in the
    // flux computation.
    if (m_proj->size() == 1)  {

        // Get pointer on HealPix projection
        const GHealpix* healpix = static_cast<const GHealpix*>(projection());

        // Get centre
        GSkyDir centre = healpix->pix2dir(pixel);

        // Get boundaries
        GSkyDirs boundaries = healpix->boundaries(pixel, 3);

        // Compute intensities
        double i0  = this->operator()(centre,         map);
        double i1  = this->operator()(boundaries[0],  map);
        double i2  = this->operator()(boundaries[4],  map);
        double i3  = this->operator()(boundaries[8],  map);
        double i4  = this->operator()(boundaries[1],  map);
        double i5  = this->operator()(boundaries[5],  map);
        double i6  = this->operator()(boundaries[9],  map);
        double i7  = this->operator()(boundaries[2],  map);
        double i8  = this->operator()(boundaries[6],  map);
        double i9  = this->operator()(boundaries[10], map);
        double i10 = this->operator()(boundaries[3],  map);
        double i11 = this->operator()(boundaries[7],  map);
        double i12 = this->operator()(boundaries[11], map);
 
        // Compute fluxes in wedges
        double flux1  = gammalib::onethird * (i0 + i1 + i2) *
                        solidangle(centre, boundaries[0], boundaries[4]);
        double flux2  = gammalib::onethird * (i0 + i2 + i3) *
                        solidangle(centre, boundaries[4], boundaries[8]);
        double flux3  = gammalib::onethird * (i0 + i3 + i4) *
                        solidangle(centre, boundaries[8], boundaries[1]);
        double flux4  = gammalib::onethird * (i0 + i4 + i5) *
                        solidangle(centre, boundaries[1], boundaries[5]);
        double flux5  = gammalib::onethird * (i0 + i5 + i6) *
                        solidangle(centre, boundaries[5], boundaries[9]);
        double flux6  = gammalib::onethird * (i0 + i6 + i7) *
                        solidangle(centre, boundaries[9], boundaries[2]);
        double flux7  = gammalib::onethird * (i0 + i7 + i8) *
                        solidangle(centre, boundaries[2], boundaries[6]);
        double flux8  = gammalib::onethird * (i0 + i8 + i9) *
                        solidangle(centre, boundaries[6], boundaries[10]);
        double flux9  = gammalib::onethird * (i0 + i9 + i10) *
                        solidangle(centre, boundaries[10], boundaries[3]);
        double flux10 = gammalib::onethird * (i0 + i10 + i11) *
                        solidangle(centre, boundaries[3], boundaries[7]);
        double flux11 = gammalib::onethird * (i0 + i11 + i12) *
                        solidangle(centre, boundaries[7], boundaries[11]);
        double flux12 = gammalib::onethird * (i0 + i12 + i1) *
                        solidangle(centre, boundaries[11], boundaries[0]);

        // Sum up fluxes
        flux = (flux1 + flux2  + flux3  + flux4 +
                flux5 + flux6  + flux7  + flux8 +
                flux9 + flux10 + flux11 + flux12);

    } // endif: we had a HealPix map

    // ... otherwise perform flux computation for WCS map. The pixel is
    // divided into 8 wedges and the flux is integrated in each wedge by
    // computing the average of the corner intensities multiplied by the
    // solid angle of the wedge. This leads to a precision of the order
    // of <1% in the flux computation.
    else {

        // Get centre
        GSkyDir centre = pix2dir(pixel);

        // Set upper pixel edge just beyond the upper boundary since the
        // upper boundary is excluded
        const double up = 0.4999999;

        // Get boundaries
        GSkyDir boundary1 = pix2dir(GSkyPixel(pixel.x()-0.5, pixel.y()-0.5));
        GSkyDir boundary2 = pix2dir(GSkyPixel(pixel.x(),     pixel.y()-0.5));
        GSkyDir boundary3 = pix2dir(GSkyPixel(pixel.x()+up,  pixel.y()-0.5));
        GSkyDir boundary4 = pix2dir(GSkyPixel(pixel.x()+up,  pixel.y()));
        GSkyDir boundary5 = pix2dir(GSkyPixel(pixel.x()+up,  pixel.y()+up));
        GSkyDir boundary6 = pix2dir(GSkyPixel(pixel.x(),     pixel.y()+up));
        GSkyDir boundary7 = pix2dir(GSkyPixel(pixel.x()-0.5, pixel.y()+up));
        GSkyDir boundary8 = pix2dir(GSkyPixel(pixel.x()-0.5, pixel.y()));

        // Compute intensities
        double i0  = this->operator()(centre,    map);
        double i1  = this->operator()(boundary1, map);
        double i2  = this->operator()(boundary2, map);
        double i3  = this->operator()(boundary3, map);
        double i4  = this->operator()(boundary4, map);
        double i5  = this->operator()(boundary5, map);
        double i6  = this->operator()(boundary6, map);
        double i7  = this->operator()(boundary7, map);
        double i8  = this->operator()(boundary8, map);

        // Compute fluxes
        double flux1 = gammalib::onethird * (i1 + i2 + i0) *
                       solidangle(boundary1, boundary2, centre);
        double flux2 = gammalib::onethird * (i2 + i3 + i0) *
                       solidangle(boundary2, boundary3, centre);
        double flux3 = gammalib::onethird * (i3 + i4 + i0) *
                       solidangle(boundary3, boundary4, centre);
        double flux4 = gammalib::onethird * (i4 + i5 + i0) *
                       solidangle(boundary4, boundary5, centre);
        double flux5 = gammalib::onethird * (i5 + i6 + i0) *
                       solidangle(boundary5, boundary6, centre);
        double flux6 = gammalib::onethird * (i6 + i7 + i0) *
                       solidangle(boundary6, boundary7, centre);
        double flux7 = gammalib::onethird * (i7 + i8 + i0) *
                       solidangle(boundary7, boundary8, centre);
        double flux8 = gammalib::onethird * (i8 + i1 + i0) *
                       solidangle(boundary8, boundary1, centre);

        // Sum up fluxes
        flux = (flux1 + flux2  + flux3  + flux4 +
                flux5 + flux6  + flux7  + flux8);

    } // endelse: we had a WCS map

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns flux within sky region
 *
 * @param[in] region Sky region.
 * @param[in] map Map index [0,...,nmaps()-1].
 * @return Flux in sky region.
 *
 * @exception GException::out_of_range
 *            Map index out of valid range
 *
 * The method returns the total flux in all sky map pixels that fall within
 * the specified region. Containment is tested using the
 * GSkyRegion::contains() method. Sky map values are multiplied by the solid
 * angle of the pixel.
 ***************************************************************************/
double GSkyMap::flux(const GSkyRegion& region, const int& map) const
{
    // Throw an error if map index is not in valid range
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_FLUX3, "Map index",
                                       map, m_num_maps);
    }

    // Initialise flux
    double flux = 0.0;

    // Loop over all sky pixels and sum flux of all pixels that fall within
    // the sky region
    for (int i = 0, ipix = m_num_pixels * map; i < m_num_pixels; ++i, ++ipix) {
        if (region.contains(inx2dir(i))) {
            flux += m_pixels(ipix) * solidangle(i);
        }
    }

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns flux within sky regions
 *
 * @param[in] regions Sky regions.
 * @param[in] map Map index [0,...,nmaps()-1].
 * @return Flux in sky region.
 *
 * @exception GException::out_of_range
 *            Map index out of valid range
 *
 * The method returns the total flux in all sky map pixels that fall within
 * the specified regions. Containment is tested using the
 * GSkyRegions::contains() method. Sky map values are multiplied by the solid
 * angle of the pixel.
 ***************************************************************************/
double GSkyMap::flux(const GSkyRegions& regions, const int& map) const
{
    // Throw an error if map index is not in valid range
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_FLUX4, "Map index",
                                       map, m_num_maps);
    }

    // Initialise flux
    double flux = 0.0;

    // Loop over all sky pixels and sum flux of all pixels that fall within
    // the sky region
    for (int i = 0, ipix = m_num_pixels * map; i < m_num_pixels; ++i, ++ipix) {
        if (regions.contains(inx2dir(i))) {
            flux += m_pixels(ipix) * solidangle(i);
        }
    }

    // Return flux
    return flux;
}


/***********************************************************************//**
 * @brief Returns array with total flux for sky maps
 *
 * @return Array of pixel flux sums.
 *
 * For a set of n sky maps, the method returns a 1D array with n entries
 * that correspond to the sum of the pixel flux in each map.
 ***************************************************************************/
GNdarray GSkyMap::flux(void) const
{
    // Initialise output array
    GNdarray flux_array = GNdarray(m_num_maps);

    // Loop over maps and store sums in output array
    for (int i = 0, ipix = 0; i < m_num_maps; ++i) {

        // Initialise flux for this map
        double flux = 0.0;

        // Compute map flux
        for (int k = 0; k < m_num_pixels; ++k, ++ipix) {
                double solidangle = this->solidangle(k);
                flux             += m_pixels(ipix) * solidangle;
        }

        // Store map flux
        flux_array(i) = flux;

    } // endfor: looped over all maps

    // Return flux array
    return flux_array;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] index Pixel index [0,...,npix()-1].
 * @return Solid angle (steradians)
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns the solid angle of the pixel with the specified @p index.
 ***************************************************************************/
double GSkyMap::solidangle(const int& index) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_SOLIDANGLE1, msg);
    }

    // Determine solid angle from pixel index.
    double solidangle = this->solidangle(inx2pix(index));

    // Return solid angle
    return solidangle;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pixel Sky map pixel.
 * @return Solid angle (steradians)
 *
 * @exception GException::invalid_value
 *            No valid sky projection found.
 *
 * Returns the solid angle of the specified sky map @p pixel.
 ***************************************************************************/
double GSkyMap::solidangle(const GSkyPixel& pixel) const
{
    // Throw error if WCS is not valid
    if (m_proj == NULL) {
        std::string msg = "Sky projection has not been defined.";
        throw GException::invalid_value(G_SOLIDANGLE2, msg);
    }

    // Initialise solid angle
    double solidangle = 0.0;

    // If pixel size matches the projection size then perform a straight
    // forward solid angle determination
    if (m_proj->size() == pixel.size()) {
        solidangle = m_proj->solidangle(pixel);
    }

    // ... otherwise, if we have a 2D projection but a 1D pixel then
    // interpret the pixel as the linear index in the pixel array
    else if (m_proj->size() == 2) {
        solidangle = m_proj->solidangle(GSkyPixel(inx2pix(int(pixel))));
    }

    // ... otherwise we have a 1D projection but a 2D pixel. There is
    // no unambiguous way to handle this case, hence we throw an exception
    else {
        std::string msg = "A 2-dimensional sky map pixel "+pixel.print()+
                          " is used to determine the solid angle for"
                          " the 1-dimensional sky projection \""+
                          m_proj->name()+"\"\n"
                          "Please specify a 1-dimensional sky map pixel.";
        throw GException::invalid_argument(G_SOLIDANGLE2, msg);
    }

    // Return solid angle
    return solidangle;
}


/***********************************************************************//**
 * @brief Returns solid angle within sky region
 *
 * @param[in] region Sky region.
 * @return Solid angle within sky region.
 *
 * The method returns the total solid angle of all sky map pixels that fall
 * within the specified region. Containment is tested using the
 * GSkyRegion::contains() method.
 ***************************************************************************/
double GSkyMap::solidangle(const GSkyRegion& region) const
{
    // Initialise flux
    double solidangle = 0.0;

    // Loop over all sky pixels and sum solid angle of all pixels that fall
    // within the sky region
    for (int i = 0; i < m_num_pixels; ++i) {
        if (region.contains(inx2dir(i))) {
            solidangle += this->solidangle(i);
        }
    }

    // Return solidangle
    return solidangle;
}


/***********************************************************************//**
 * @brief Returns solid angle within sky regions
 *
 * @param[in] regions Sky regions.
 * @return Solid angle within sky region.
 *
 * The method returns the total solid angle of all sky map pixels that fall
 * within the specified regions. Containment is tested using the
 * GSkyRegions::contains() method.
 ***************************************************************************/
double GSkyMap::solidangle(const GSkyRegions& regions) const
{
    // Initialise flux
    double solidangle = 0.0;

    // Loop over all sky pixels and sum solid angle of all pixels that fall
    // within the sky regions
    for (int i = 0; i < m_num_pixels; ++i) {
        if (regions.contains(inx2dir(i))) {
            solidangle += this->solidangle(i);
        }
    }

    // Return solidangle
    return solidangle;
}


/***********************************************************************//**
 * @brief Set sky projection
 *
 * @param[in] proj Sky projection.
 *
 * Sets the projection from celestial to pixel coordinates. The method
 * performs a deep copy of @p proj, allowing to destroy the argument after
 * using the method.
 *
 * Warning: this method may corrupt the GSkyMap object as it allows assigning
 * for example a 1D projection to a 2D skymap. Please use this method only
 * when you know what you're doing.
 *
 * @todo We may restrict this method to not allow changing the projection
 * dimension.
 ***************************************************************************/
void GSkyMap::projection(const GSkyProjection& proj)
{
    // Free existing projection only if it differs from current projection.
    // This prevents unintential deallocation of the argument
    if ((m_proj != NULL) && (m_proj != &proj)) {
        delete m_proj;
    }

    // Clone input WCS
    m_proj = proj.clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks if sky direction falls in map
 *
 * @param[in] dir Sky direction.
 * @return True if sky direction falls into map, false otherwise.
 *
 * This method checks if the specified sky direction falls within the pixels
 * covered by the skymap. The method uses the dir2xy method to convert the
 * sky direction into 2D pixel indices, and then checks whether the pixel
 * indices fall in the skymap.
 ***************************************************************************/
bool GSkyMap::contains(const GSkyDir& dir) const
{
    // Initialise containment
    bool contains = false;

    // Convert sky direction into sky pixel
    try {

        // Convert sky direction into sky pixel
        GSkyPixel pixel = dir2pix(dir);

        // Check pixel containment
        contains = this->contains(pixel);

    }
    catch (GException::invalid_argument) {
        contains = false;
    }

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Checks if sky map pixel falls in map
 *
 * @param[in] pixel Sky map pixel.
 * @return Trus if pixels is within map, false otherwise.
 *
 * Checks whether the specified sky map @p pixel falls within the skymap
 * or not. A pixel is considered to fall in the sky may if the value is
 * contained in the interval
 *
 * \f[
 *    [-0.5, n-0.5[
 * \f]
 *
 * where \f$n\f$ is the number of sky map pixels. Note that the upper bound
 * is excluded from the interval.
 ***************************************************************************/
bool GSkyMap::contains(const GSkyPixel& pixel) const
{
    // Initialise containment flag
    bool inmap = false;

    // Test 1D pixels
    if (pixel.is_1D()) {
        if (pixel.index() >= -0.5 &&
            pixel.index() <  double(m_num_pixels)-0.5) {
            inmap = true;
        }
    }

    // Test 2D pixels
    else if (pixel.is_2D()) {

        // If pixel is in range then set containment flag to true
        if ((pixel.x() >= -0.5 && pixel.x() < double(m_num_x)-0.5) &&
            (pixel.y() >= -0.5 && pixel.y() < double(m_num_y)-0.5)) {
            inmap = true;
        }

    }

    // Return containment flag
    return inmap;
}


/***********************************************************************//**
 * @brief Checks whether a region overlaps with this map
 *
 * @param[in] region Region
 * @return True if region overlaps with map, false otherwise
 *
 * @exception GException::feature_not_implemented
 *            Region is not a circular sky region
 ***************************************************************************/
bool GSkyMap::overlaps(const GSkyRegion& region) const
{
    // Initialise overlap
    bool overlap = false;

    // If the region is a circular region than test on overlap
    if (region.type() == "Circle") {
        const GSkyRegionCircle* circle =
              static_cast<const GSkyRegionCircle*>(&region);
        overlap = overlaps_circle(*circle);
    }

    // ... otherwise signal that method is not implemented
    else {
        std::string msg = "Method not implemented for non-circular regions.";
        throw GException::feature_not_implemented(G_OVERLAPS, msg);
    }

    // Return overlap flag
    return overlap;
}


/***********************************************************************//**
 * @brief Extract maps into a new sky map object
 *
 * @param[in] map First map to extract
 * @param[in] nmaps Number of maps to extract
 * @return Extracted map(s).
 *
 * @exception GException::out_of_range
 *            First map index outside valid range
 * @exception GException::invalid_argument
 *            Requested number of maps are not available.
 *
 * Extracts @p nmaps sky maps starting from @p map from the sky map.
 ***************************************************************************/
GSkyMap GSkyMap::extract(const int& map, const int& nmaps) const
{
    // Throw an exception if the first map index is invalid
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_EXTRACT, "Sky map index", map,
                                       m_num_maps);
    }

    // Throw an exception if the number of maps is invalid
    if (nmaps < 0) {
        std::string msg = "The number of maps to extract cannot be negative "
                          "(nmaps="+gammalib::str(nmaps)+").";
        throw GException::invalid_argument(G_EXTRACT, msg);
    }
    if (nmaps > m_num_maps-map) {
        std::string msg = "The number of maps to extract ("+
                          gammalib::str(nmaps)+") exceeds the number of maps "
                          "that are available ("+
                          gammalib::str(m_num_maps-map)+").";
        throw GException::invalid_argument(G_EXTRACT, msg);
    }

    // Allocate extracted pixels
    GNdarray extracted_pixels(m_num_pixels, nmaps);

    // Extract pixels
    int           n_size = m_num_pixels * nmaps;
    const double *src    = m_pixels.data() + map*m_num_pixels;
    double       *dst    = extracted_pixels.data();
    for (int i = 0; i < n_size; ++i) {
        *dst++ = *src++;
    }

    // Create a copy of the map
    GSkyMap result = *this;

    // Attach copied pixels to the map
    result.m_pixels = extracted_pixels;

    // Set number of maps
    result.m_num_maps = nmaps;

    // Reset shape of maps since we do not know how to arrange the maps
    result.m_shape.clear();
    result.m_shape.push_back(nmaps);

    // Return map
    return result;
}


/***********************************************************************//**
 * @brief Extract a sub-range of pixels in the map in x,y
 *
 * @param[in] startx    Starting bin in X (inclusive)
 * @param[in] stopx     Last bin in X (inclusive)
 * @param[in] starty    Starting bin in Y (inclusive)
 * @param[in] stopy     Last bin in Y (inclusive)
 * @return Skymap that overlaps with supplied exclusions
 *
 * @exception GException::invalid_argument
 *            Method not valid for HPX projection
 *
 * This method creates a new skymap consisting of all pixels in the map in the 
 * range [startx,stopx] and [starty,stopy]. The boundary values provided are 
 * inclusive and the actual projection is unchanged. The number of maps in the 
 * skymap is also unchanged.
 * 
 * NOTE: The pixel ranges are indexed from 0, as in
 ***************************************************************************/
GSkyMap GSkyMap::extract(const int& startx, const int& stopx,
                         const int& starty, const int& stopy) const
{
    // Make sure this isn't a 'HealPix' map
    if (m_proj->code() == "HPX") {
        std::string msg = "Method not valid for \"HPX\" projection. Please "
                          "specify WCS projection.";
        throw GException::invalid_argument(G_EXTRACT_INT, msg);
    }
    else if ((startx > stopx) || (starty > stopy)) {
        throw GException::invalid_argument(G_EXTRACT_INT,
                                           "Invalid x,y range provided");
    }

    // Define the actual range of pixels
    int first_x = (startx < 0)       ? 0         : startx;
    int last_x  = (stopx >= m_num_x) ? m_num_x-1 : stopx;
    int first_y = (starty < 0)       ? 0         : starty;
    int last_y  = (stopy >= m_num_y) ? m_num_y-1 : stopy;
    int nxpix   = last_x - first_x + 1;
    int nypix   = last_y - first_y + 1;
    int npixels = nxpix * nypix;    

    // Create an array of values to be plotted
    GNdarray new_pixels(nxpix, nypix, m_num_maps);

    // Get new pixel values
    for (int ix = 0; ix < nxpix; ++ix) {
        int ix_src = ix + first_x;
        for (int iy = 0; iy < nypix; ++iy) {
            int iy_src = iy + first_y;
            for (int imap = 0; imap < m_num_maps; ++imap) {
                new_pixels(ix, iy, imap) = m_pixels(ix_src, iy_src, imap);
            }
        }
    }

    // Create a copy of this map
    GSkyMap result = *this;

    // Update the projection
    GWcs* proj = static_cast<GWcs*>(result.m_proj);
    proj->set(proj->coordsys(),
              proj->crval(0), proj->crval(1),
              proj->crpix(0) - first_x,             // Offset reference pixel
              proj->crpix(1) - first_y,             // Offset reference pixel
              proj->cdelt(0), proj->cdelt(1));

    // Update the pixel information
    result.m_pixels     = new_pixels;
    result.m_num_x      = nxpix;
    result.m_num_y      = nypix;
    result.m_num_pixels = npixels;

    // Re-initialise computation cache
    result.m_hascache  = false;
    result.m_contained = false;
    result.m_last_dir.clear();
    result.m_interpol.clear();

    // Return map
    return result;
}


/***********************************************************************//**
 * @brief Extract the spatial portion of the maps that overlap @p inclusions
 *
 * @param[in] inclusions List of GSkyRegion objects
 * @return Skymap that overlaps with supplied exclusions
 *
 * @exception GException::invalid_argument
 *            Method not valid for HPX projection
 *
 * This method computes the x,y range that covers the regions in @p inclusions, 
 * then returns that range as a new GSkyMap.
 * 
 * NOTE: A 1-pixel border is added to the returned map to provide better
 * coverage in the returned map.
 ***************************************************************************/
GSkyMap GSkyMap::extract(const GSkyRegions& inclusions) const
{
    // Make sure this isn't a 'HealPix' map
    if (m_proj->code() == "HPX") {
        std::string msg = "Method not valid for \"HPX\" projection. Please "
                          "specify WCS projection.";
        throw GException::invalid_argument(G_EXTRACT_REG, msg);
    }

    // Preliminary variables
    int startx = nx()+1;
    int stopx  = 0;
    int starty = ny()+1;
    int stopy  = 0;

    // Loop on pixels
    for (int i = 0; i < npix(); ++i) {

        // Get the pixel and direction
        GSkyPixel pixel  = inx2pix(i);
        GSkyDir   pixdir = pix2dir(pixel);

        // Check if the pixel is covered by 'inclusions'
        if (inclusions.contains(pixdir)) {

            // Get the x,y pixel
            int pix_x = int(pixel.x());
            int pix_y = int(pixel.y());

            // Update x,y pixel ranges
            startx = (pix_x < startx) ? pix_x : startx;
            stopx  = (pix_x > stopx)  ? pix_x : stopx;
            starty = (pix_y < starty) ? pix_y : starty;
            stopy  = (pix_y > stopy)  ? pix_y : stopy;

        } // endif: pixel was covered by inclusions

    } // endfor: looped over all pixels

    // Trim the skymap with 1 pixel buffer all around.
    return (this->extract(startx-1, stopx+1, starty-1, stopy+1));
}


/***********************************************************************//**
 * @brief Stack all maps into a single map
 *
 * The methods replaces the sky map by a version with only a single map
 * by summing over the pixel values for all maps. If the sky map has no
 * pixels or there is only a single map in the object, the method does
 * nothing.
 ***************************************************************************/
void GSkyMap::stack_maps(void)
{
    // Continue only if the map has pixels and if there is more than one map
    if (m_num_pixels > 0 && m_num_maps > 1) {

        // Initialise stacked pixels
        GNdarray stacked_pixels;

        // Allocate stacked pixels
        if (m_proj->code() == "HPX") {
            stacked_pixels = GNdarray(m_num_pixels);
        }
        else {
            stacked_pixels = GNdarray(m_num_x, m_num_y);
        }

        // Stack map and save in pixels
        for (int i = 0; i < m_num_pixels; ++i) {
            double sum = 0.0;
            for (int imap = 0; imap < m_num_maps; ++imap) {
                sum += (*this)(i,imap);
            }
            stacked_pixels(i) = sum;
        }

        // Set pixels to stacked pixels
        m_pixels = stacked_pixels;

        // Set number of maps to 1
        m_num_maps = 1;

        // Reset shape of maps
        m_shape.clear();
        m_shape.push_back(m_num_maps);

    } // endif: map had pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return sky region circle that encloses the sky map
 *
 * @return Enclosing sky region circle.
 *
 * Returns a sky region circle that encloses the sky map.
 *
 * For a sky map in HealPix projection the method returns a region circle
 * that encloses the full sky, with a centre at Right Ascension and
 * Declination of 0 and a radius of 180 degrees.
 *
 * For a sky map in World Coordinate projection the method returns a region
 * circle that is centred on the central pixel of the sky map (i.e. the
 * pixel with indices \f$(n_x/2,n_y/2)\f$, where \f$n_x\f$ and \f$n_y\f$ are
 * the number of pixels in the X- and Y-drection) and a radius that is
 * determined from the maximum radial distance of all pixels with respect
 * to the central pixel.
 ***************************************************************************/
GSkyRegionCircle GSkyMap::region_circle(void) const
{
    // Initialise sky region circle
    GSkyRegionCircle region;

    // If the sky map is in HealPix projection then set radius to 180 deg
    if (projection()->code() == "HPX") {
        region = GSkyRegionCircle(0.0, 0.0, 180.0);
    }

    // ... otherwise sky map is in World Coordinate Projection
    else {

        // Initialise maximum radius
        double max_radius = 0.0;
        int    max_index  =  -1;

        // Determine map centre
        GSkyPixel pixel(nx()/2.0, ny()/2.0);
        GSkyDir   centre = pix2dir(pixel);

        // Determine maximum pixel distance with respect to map centre
        for (int i = 0; i < npix(); ++i) {
            double radius = inx2dir(i).dist_deg(centre);
            if (radius > max_radius) {
                max_radius = radius;
                max_index  = i;
            }
        }

        // Consider edges of pixel with maximum distance
        if (max_index != -1) {

            // Set upper pixel edge just beyond the upper boundary since the
            // upper boundary is excluded
            const double up = 0.4999999;

            // Get pixel with maximum distance
            GSkyPixel pixel = inx2pix(max_index);

            // Get boundaries of pixel
            GSkyDir boundary1 = pix2dir(GSkyPixel(pixel.x()-0.5, pixel.y()-0.5));
            GSkyDir boundary2 = pix2dir(GSkyPixel(pixel.x(),     pixel.y()-0.5));
            GSkyDir boundary3 = pix2dir(GSkyPixel(pixel.x()+up,  pixel.y()-0.5));
            GSkyDir boundary4 = pix2dir(GSkyPixel(pixel.x()+up,  pixel.y()));
            GSkyDir boundary5 = pix2dir(GSkyPixel(pixel.x()+up,  pixel.y()+up));
            GSkyDir boundary6 = pix2dir(GSkyPixel(pixel.x(),     pixel.y()+up));
            GSkyDir boundary7 = pix2dir(GSkyPixel(pixel.x()-0.5, pixel.y()+up));
            GSkyDir boundary8 = pix2dir(GSkyPixel(pixel.x()-0.5, pixel.y()));

            // Get radii
            double radius1 = boundary1.dist_deg(centre);
            double radius2 = boundary2.dist_deg(centre);
            double radius3 = boundary3.dist_deg(centre);
            double radius4 = boundary4.dist_deg(centre);
            double radius5 = boundary5.dist_deg(centre);
            double radius6 = boundary6.dist_deg(centre);
            double radius7 = boundary7.dist_deg(centre);
            double radius8 = boundary8.dist_deg(centre);

            // Determine maximum radius
            if (radius1 > max_radius) {
                max_radius = radius1;
            }
            if (radius2 > max_radius) {
                max_radius = radius2;
            }
            if (radius3 > max_radius) {
                max_radius = radius3;
            }
            if (radius4 > max_radius) {
                max_radius = radius4;
            }
            if (radius5 > max_radius) {
                max_radius = radius5;
            }
            if (radius6 > max_radius) {
                max_radius = radius6;
            }
            if (radius7 > max_radius) {
                max_radius = radius7;
            }
            if (radius8 > max_radius) {
                max_radius = radius8;
            }

        } // endif: considered edges of pixel

        // Set sky region
        region = GSkyRegionCircle(centre, max_radius);

    } // endelse: sky map was in World Coordinate Projection

    // Return region
    return region;
}


/***********************************************************************//**
 * @brief Load skymap from FITS file.
 *
 * @param[in] filename FITS file name..
 *
 * Loads HEALPix and non HEALPix skymaps. First searches for HEALPix map in
 * FITS file by scanning all HDUs for PIXTYPE=HEALPIX. If no HEALPix map has
 * been found then search load first non-empty image.
 *
 * @todo Do we have to restrict a HEALPix map to a BinTable and a WCS map
 * to a Double precision image???
 ***************************************************************************/
void GSkyMap::load(const GFilename& filename)
{
    // Clear sky map
    clear();

    // Initialize load flag
    bool loaded(false);

    // Open FITS file
    GFits fits(filename);

    // If an extension name is specified then first try loading that
    // extension
    if (filename.has_extname()) {
        const GFitsHDU& hdu = *fits.at(filename.extname());
        if (is_healpix(hdu)) {
            read_healpix(static_cast<const GFitsTable&>(hdu));
            loaded = true;
        }
        else if (is_wcs(hdu)) {
            read_wcs(static_cast<const GFitsImage&>(hdu));
            loaded = true;
        }
    }

    // ... otherwise is an extension number is specified then try loading
    // the corresponding HDU
    else if (filename.has_extno()) {
        const GFitsHDU& hdu = *fits.at(filename.extno());
        if (is_healpix(hdu)) {
            read_healpix(static_cast<const GFitsTable&>(hdu));
            loaded = true;
        }
        else if (is_wcs(hdu)) {
            read_wcs(static_cast<const GFitsImage&>(hdu));
            loaded = true;
        }
    }

    // If no map has yet been loaded then scan the file for an appropriate
    // HDU and read sky map from the first suitable HDU
    if (!loaded) {

        // Get number of HDUs
        int num = fits.size();

        // First search for HEALPix extension. We can skip the first
        // extension since this is always an image and a HEALPix map
        // is stored in a binary table
        for (int extno = 1; extno < num; ++extno) {

            // Get reference to HDU
            const GFitsHDU& hdu = *fits.at(extno);

            // If HDU is HEALPix then read data
            if (is_healpix(hdu)) {
                read_healpix(static_cast<const GFitsTable&>(hdu));
                loaded = true;
                break;
            }

        } // endfor: looped over HDUs

        // If we have not found a HEALPIX map then search now for an
        // image.
        if (!loaded) {
            for (int extno = 0; extno < num; ++extno) {

                // Get referene to HDU
                const GFitsHDU& hdu = *fits.at(extno);

                // If HDU is WCS then read data
                if (is_wcs(hdu)) {
                    read_wcs(static_cast<const GFitsImage&>(hdu));
                    //loaded = true;
                    break;
                }

            } // endfor: looped over HDUs
        } // endif: no sky map yet loaded
    } // endif: no sky map yet loaded

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save sky map into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves the sky map into a FITS file. If the file exists already and the
 * @p clobber parameter is true, the method will overwrite the content of
 * the existing file. Otherwise, an exception is thrown.
 ***************************************************************************/
void GSkyMap::save(const GFilename& filename, const bool& clobber) const
{
    // Continue only if we have data to save
    if (m_proj != NULL) {

        // Initialise HDU pointer
        GFitsHDU* hdu = NULL;

        // Case A: Skymap is Healpix
        if (m_proj->code() == "HPX") {
            hdu = create_healpix_hdu();
        }

        // Case B: Skymap is not Healpix
        else {
            hdu = create_wcs_hdu();
        }

        // If we have a valid HDU then save it now to the FITS file
        if (hdu != NULL) {

            // Set extension name
            if (filename.has_extname()) {
                hdu->extname(filename.extname());
            }

            // Initialise empty FITS file
            GFits fits;

            // Append sky map to FITS file
            fits.append(*hdu);

            // Save FITS file (without extension name which was extracted
            // earlier and set in the HDU)
            fits.saveto(filename.url(), clobber);

            // Delete HDU
            delete hdu;

        } // endif: HDU was valid

    } // endif: we had data to save

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read skymap from FITS HDU
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
void GSkyMap::read(const GFitsHDU& hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // If HDU is a Healpix HDU then load a Healpix map
    if (is_healpix(hdu)) {
        read_healpix(static_cast<const GFitsTable&>(hdu));
    }

    // ... otherwise try loading as WCS map
    else if (is_wcs(hdu)) {
        read_wcs(static_cast<const GFitsImage&>(hdu));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write sky map into FITS file
 *
 * @param[in] file FITS file pointer.
 * @param[in] extname Sky map extension name.
 * @return Pointer to written HDU.
 *
 * Write the sky map into a FITS file. Optionally, the extension name of
 * the FITS HDU can be specified using the @p extname parameter. The method
 * returns a pointer to the appended HDU. If no HDU has been appended, for
 * example because the sky map is empty, the method will return NULL.
 ***************************************************************************/
GFitsHDU* GSkyMap::write(GFits& file, const std::string& extname) const
{
    // Initialise pointer to appended HDU
    GFitsHDU* hdu_appended = NULL;

    // Continue only if we have data to save
    if (m_proj != NULL) {

        // Initialise local HDU pointer
        GFitsHDU* hdu = NULL;

        // Case A: Skymap is Healpix
        if (m_proj->code() == "HPX") {
            hdu = create_healpix_hdu();
        }

        // Case B: Skymap is not Healpix
        else {
            hdu = create_wcs_hdu();
        }

        // Append HDU to FITS file.
        if (hdu != NULL) {

            // Optionally set extension name
            if (extname.length() > 0) {
                hdu->extname(extname);
            }

            // Append HDU and recover the pointer to the appended HDU
            hdu_appended = file.append(*hdu);

            // Delete HDU
            delete hdu;

        } // endif: there was a valid HDU

    } // endif: we had data to save

    // Return pointer to appended HDU
    return hdu_appended;
}


/***********************************************************************//**
 * @brief Publish sky map
 *
 * @param[in] name Name of sky map.
 *
 * Publishes the sky map on a Virtual Observatory Hub. If no Hub is currently
 * active, the method will start a new Hub. If on sky map has been allocated
 * the method does nothing.
 ***************************************************************************/
void GSkyMap::publish(const std::string& name) const
{
    // Create FITS file containing the sky map
    GFits fits;

    // Write sky map into FITS file and return a pointer to the written HDU.
    // This method returns a NULL pointer if no sky map projection exists,
    // which typically is true for empty sky maps.
    GFitsHDU* hdu = write(fits);

    // Continue only if HDU pointer is valid
    if (hdu != NULL) {

        // Optionally set extension name
        if (!name.empty()) {
            hdu->extname(name);
        }

        // Create VO Client
        GVOClient client;

        // Publish map using VO client
        client.publish(*hdu);

    } // endif: HDU pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print sky map
 *
 * @param[in] chatter Chattiness.
 * @return String containing sky map information.
 ***************************************************************************/
std::string GSkyMap::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSkyMap ===");

        // Append number of pixels and number of maps
        result.append("\n"+gammalib::parformat("Number of pixels"));
        result.append(gammalib::str(m_num_pixels));

        // Append WCS dimension information
        if (m_proj != NULL && m_proj->size() == 2) {
            result.append("\n"+gammalib::parformat("X axis dimension"));
            result.append(gammalib::str(m_num_x));
            result.append("\n"+gammalib::parformat("Y axis dimension"));
            result.append(gammalib::str(m_num_y));
        }

        // Append number of maps
        result.append("\n"+gammalib::parformat("Number of maps"));
        result.append(gammalib::str(m_num_maps));

        // Append shape of maps
        result.append("\n"+gammalib::parformat("Shape of maps"));
        if (m_shape.size() > 0) {
            std::string shape = "(";
            for (int i = 0; i < m_shape.size(); ++i) {
                if (i > 0) {
                    shape += ", ";
                }
                shape += gammalib::str(m_shape[i]);
            }
            shape += ")";
            result.append(shape);
        }
        else {
            result.append("not defined");
        }

        // Append sky projection information
        if (gammalib::reduce(chatter) > SILENT) {
            if (m_proj != NULL) {
                result.append("\n"+m_proj->print(gammalib::reduce(chatter)));
            }
            else {
                result.append("\n"+gammalib::parformat("Sky projection"));
                result.append("not defined");
            }
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyMap::init_members(void)
{
    // Initialise members
    m_num_pixels = 0;
    m_num_maps   = 0;
    m_num_x      = 0;
    m_num_y      = 0;
    m_shape.clear();
    m_proj       = NULL;
    m_pixels.clear();

    // Initialise computation cache
    m_hascache    = false;
    m_contained   = false;
    m_last_dir.clear();
    m_interpol.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] map Sky map.
 ***************************************************************************/
void GSkyMap::copy_members(const GSkyMap& map)
{
    // Copy attributes
    m_num_pixels = map.m_num_pixels;
    m_num_maps   = map.m_num_maps;
    m_num_x      = map.m_num_x;
    m_num_y      = map.m_num_y;
    m_shape      = map.m_shape;
    m_pixels     = map.m_pixels;

    // Copy computation cache
    m_hascache  = map.m_hascache;
    m_contained = map.m_contained;
    m_last_dir  = map.m_last_dir;
    m_interpol  = map.m_interpol;

    // Clone sky projection if it is valid
    if (map.m_proj != NULL) m_proj = map.m_proj->clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyMap::free_members(void)
{
    // Free memory
    if (m_proj != NULL) {
        delete m_proj;
    }

    // Signal free pointers
    m_proj = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set World Coordinate System
 *
 * @param[in] wcs World Coordinate System code.
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel.
 * @param[in] crpix2 Y index of reference pixel.
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 *
 * @exception GException::invalid_argument
 *            Invalid wcs parameter (World Coordinate System not known).
 *
 * This method sets the WCS projection pointer based on the WCS code and
 * sky map parameters. It makes use of the GWcsRegistry class to allocate
 * the correct derived class. Note that this method does not support the
 * HPX projection.
 ***************************************************************************/
void GSkyMap::set_wcs(const std::string& wcs,
                      const std::string& coords,
                      const double&      crval1,
                      const double&      crval2,
                      const double&      crpix1,
                      const double&      crpix2,
                      const double&      cdelt1,
                      const double&      cdelt2)
{
    // Convert WCS to upper case
    std::string uwcs = gammalib::toupper(wcs);

    // Check if HPX was requested (since this is not allowed)
    if (uwcs == "HPX") {
        std::string msg = "Method not valid for \"HPX\" projection. Please "
                          "specify WCS projections.";
        throw GException::invalid_argument(G_SET_WCS, msg);
    }

    // ... otherwise get projection from registry
    else {
        // Allocate WCS registry
        GWcsRegistry registry;

        // Allocate projection from registry
        GWcs* projection = registry.alloc(uwcs);
        m_proj           = projection;

        // Signal if projection type is not known
        if (projection == NULL) {
            std::string msg = "WCS projection code \""+uwcs+"\" not known. "
                              "Please specify one of the following WCS "
                              "projections: "+registry.content()+".";
            throw GException::invalid_argument(G_SET_WCS, msg);
        }

        // Setup WCS
        projection->set(coords, crval1, crval2, crpix1, crpix2,
                        cdelt1, cdelt2);

    } // endelse: got projection from registry

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Healpix data from FITS table.
 *
 * @param[in] table FITS table.
 *
 * @exception GException::runtime_error
 *            Incompatible HealPix projection encountered.
 *
 * HEALPix data may be stored in various formats depending on the 
 * application that has writted the data. HEALPix IDL, for example, may
 * store the data in vectors of length 1024 if the number of pixels is
 * a multiple of 1024. On the other hand, vectors may also be used to store
 * several HEALPix maps into a single column. Alternatively, multiple maps
 * may be stored in multiple columns.
 ***************************************************************************/
void GSkyMap::read_healpix(const GFitsTable& table)
{
    // Determine number of rows and columns in table
    int nrows = table.nrows();
    int ncols = table.ncols();
    #if defined(G_READ_HEALPIX_DEBUG)
    std::cout << "nrows=" << nrows << " ncols=" << ncols << std::endl;
    #endif

    // Allocate Healpix projection
    m_proj = new GHealpix;

    // Read WCS information from FITS header
    m_proj->read(table);

    // Set number of pixels based on NSIDE parameter
    m_num_pixels = static_cast<GHealpix*>(m_proj)->npix();
    #if defined(G_READ_HEALPIX_DEBUG)
    std::cout << "m_num_pixels=" << m_num_pixels << std::endl;
    #endif

    // Number of map pixels has to be a multiple of the number of
    // rows in column
    if (m_num_pixels % nrows != 0) {
        std::string msg = "Number of "+gammalib::str(nrows)+" rows in FITS "
                          "table must correspond to number of "+
                          gammalib::str(m_num_pixels)+" pixels in HealPix "
                          "projection.";
        throw GException::runtime_error(G_READ_HEALPIX, msg);
    }

    // Determine vector length for HEALPix data storage
    int nentry = m_num_pixels / nrows;
    #if defined(G_READ_HEALPIX_DEBUG)
    std::cout << "nentry=" << nentry << std::endl;
    #endif

    // Determine number of maps from the number of maps that fit into
    // all columns. Only count columns that can fully hold the map.
    m_num_maps = 0;
    for (int icol = 0; icol < ncols; ++icol) {
        const GFitsTableCol* col = table[icol];
        if (col->number() % nentry == 0) {
            m_num_maps += col->number() / nentry;
        }
    }
    #if defined(G_READ_HEALPIX_DEBUG)
    std::cout << "m_num_maps=" << m_num_maps << std::endl;
    #endif

    // Initialise shape of maps
    m_shape.clear();
    m_shape.push_back(m_num_maps);

    // Allocate pixels to hold the map
    m_pixels = GNdarray(m_num_pixels, m_num_maps);

    // Initialise map counter
    int imap = 0;

    // Loop over all columns
    for (int icol = 0; icol < ncols; ++icol) {

        // Get next column
        const GFitsTableCol* col = table[icol];

        // Only consider columns that can fully hold maps
        if (col->number() % nentry == 0) {

            // Determine number of maps in column
            int num = col->number() / nentry;

            // Loop over all maps in column
            int inx_start = 0;
            int inx_end   = nentry;
            for (int i = 0; i < num; ++i) {

                // Load map
                double *ptr = m_pixels.data() + m_num_pixels*imap;
                for (int row = 0; row < col->nrows(); ++row) {
                    for (int inx = inx_start; inx < inx_end; ++inx) {
                        *ptr++ = col->real(row,inx);
                    }
                }
                #if defined(G_READ_HEALPIX_DEBUG)
                std::cout << "Load map=" << imap << " index="
                          << inx_start << "-" << inx_end << std::endl;
                #endif

                // Increment index range
                inx_start  = inx_end;
                inx_end   += nentry;

                // Increment map counter
                imap++;

                // Break if we have loaded all maps
                if (imap >= m_num_maps) {
                    break;
                }

            } // endfor: looped over all maps in column
        } // endif: column could fully hold maps

        // Break if we have loaded all maps
        if (imap >= m_num_maps) {
            break;
        }

    } // endfor: looped over all columns

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read WCS image from FITS HDU
 *
 * @param[in] image FITS image.
 *
 * @exception GException::invalid_argument
 *            FITS image has less than two dimensions
 * @exception GException::invalid_value
 *            Sky map covers more than 360 deg in longitude or 180 deg in
 *            latitude
 *
 * Reads sky maps from a FITS image extension containing a set of maps given
 * in World Coordinate System. The method handles general n-dimensional
 * images and sets the map shape attribute according to the number of map
 * dimensions found in the FITS HDU.
 ***************************************************************************/
void GSkyMap::read_wcs(const GFitsImage& image)
{
    // Throw an exception if the FITS image is not at least a 2D image
    if (image.naxis() < 2) {
        std::string msg = "Sky map has "+gammalib::str(image.naxis())+
                          " dimensions, which is less than the two dimensions"
                          " that are required for a WCS image.";
        throw GException::invalid_argument(G_READ_WCS, msg);
    }

    // Allocate WCS
    alloc_wcs(image);

    // Read projection information from FITS header
    m_proj->read(image);

    // Clear map shape
    m_shape.clear();

    // Extract map dimension, number of maps and map shape from FITS image
    m_num_x = image.naxes(0);
    m_num_y = image.naxes(1);
    if (image.naxis() == 2) {
        m_num_maps = 1;
        m_shape.push_back(m_num_maps);
    }
    else {
        m_num_maps = image.naxes(2);
        m_shape.push_back(m_num_maps);
        for (int i = 3; i < image.naxis(); ++i) {
            m_num_maps *= image.naxes(i);
            m_shape.push_back(image.naxes(i));
        }
    }
    #if defined(G_READ_WCS_DEBUG)
    std::cout << "m_num_x=" << m_num_x << std::endl;
    std::cout << "m_num_y=" << m_num_y << std::endl;
    std::cout << "m_num_maps=" << m_num_maps << std::endl;
    #endif

    // Compute number of pixels
    m_num_pixels = m_num_x * m_num_y;
    #if defined(G_READ_WCS_DEBUG)
    std::cout << "m_num_pixels=" << m_num_pixels << std::endl;
    #endif

    // Allocate pixels to hold the map
    m_pixels = GNdarray(m_num_x, m_num_y, m_num_maps);

    // Read image
    double* ptr  = m_pixels.data();
    int     size = m_num_pixels * m_num_maps;
    for (int i = 0; i < size; ++i) {
        *ptr++ = (double)image.pixel(i);
    }

    // Check if the map is too large
    double x_size = std::abs(m_num_x * static_cast<GWcs*>(m_proj)->cdelt(0));
    double y_size = std::abs(m_num_y * static_cast<GWcs*>(m_proj)->cdelt(1));
    if (x_size > 360.001) {
        std::string msg = "Skymap covers "+gammalib::str(x_size)+" degrees "
                          "in the X axis which results in ambiguous "
                          "coordinate transformations. Please provide a "
                          "skymap that does not cover more than 360 degrees.";
        throw GException::invalid_value(G_READ_WCS, msg);
    }
    if (y_size > 180.001) {
        std::string msg = "Skymap covers "+gammalib::str(y_size)+" degrees "
                          "in the Y axis which results in ambiguous "
                          "coordinate transformations. Please provide a "
                          "skymap that does not cover more than 180 degrees.";
        throw GException::invalid_value(G_READ_WCS, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate WCS class
 *
 * @param[in] image FITS image.
 *
 * @exception GException::invalid_argument
 *            CTYPE1 and CTYPE2 keywords are incompatible or WCS projection
 *            of FITS file not supported by GammaLib.
 ***************************************************************************/
void GSkyMap::alloc_wcs(const GFitsImage& image)
{
    // Get standard keywords
    std::string ctype1 = image.string("CTYPE1");
    std::string ctype2 = image.string("CTYPE2");

    // Extract projection type
    std::string xproj = ctype1.substr(5,3);
    std::string yproj = ctype2.substr(5,3);

    // Check that projection type is identical on both axes
    if (xproj != yproj) {
        std::string msg = "X axis projection \""+xproj+"\" differs from Y "
                          "axis projection \""+yproj+"\". Please specify a "
                          "FITS image with a valid WCS header.";
        throw GException::invalid_argument(G_ALLOC_WCS, msg);
    }

    // Allocate WCS registry
    GWcsRegistry registry;

    // Allocate projection from registry
    m_proj = registry.alloc(xproj);

    // Signal if projection type is not known
    if (m_proj == NULL) {
        std::string msg = "WCS projection code \""+xproj+"\" not known. Please "
                          "specify a FITS image with one of the following "
                          "WCS projections: "+registry.content()+".";
        throw GException::invalid_argument(G_ALLOC_WCS, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create FITS HDU containing Healpix data
 *
 * This method allocates a binary table HDU that contains the Healpix data.
 * Deallocation of the table has to be done by the client.
 ***************************************************************************/
GFitsBinTable* GSkyMap::create_healpix_hdu(void) const
{
    // Initialise result to NULL pointer
    GFitsBinTable* hdu = NULL;

    // Compute size of Healpix data
    int size = m_num_pixels * m_num_maps;

    // Continue only if we have pixels
    if (size > 0) {

        // Set number of rows and columns
        int rows   = m_num_pixels;
        int number = m_num_maps;

        // Create column to hold Healpix data
        GFitsTableDoubleCol column = GFitsTableDoubleCol("DATA", rows, number);

        // Fill data into column
        const double* ptr = m_pixels.data();
        for (int inx = 0; inx < number; ++inx) {
            for (int row = 0; row < rows; ++row) {
                column(row,inx) = *ptr++;
            }
        }

        // Create HDU that contains Healpix map in a binary table
        hdu = new GFitsBinTable(rows);
        hdu->append(column);

    } // endif: there were pixels

    // ... otherwise create an empty header
    else {
        hdu = new GFitsBinTable;
    }

    // Set extension name
    hdu->extname("HEALPIX");

    // If we have WCS information then write into FITS header
    if (m_proj != NULL) {
        m_proj->write(*hdu);
    }

    // Set additional keywords
    hdu->card("NBRBINS", m_num_maps, "Number of HEALPix maps");

    // Return HDU
    return hdu;
}


/***********************************************************************//**
 * @brief Create FITS HDU containing WCS image
 *
 * This method allocates an image HDU that contains the WCS image data.
 * Deallocation of the image has to be done by the client.
 *
 * @todo Set additional keywords.
 ***************************************************************************/
GFitsImageDouble* GSkyMap::create_wcs_hdu(void) const
{
    // Initialise result to NULL pointer
    GFitsImageDouble* hdu = NULL;

    // Compute size of Healpix data
    int size = m_num_pixels * m_num_maps;

    // Continue only if we have pixels
    if (size > 0) {

        // Set dimension vector of all axes. In case that only one map
        // exists then create simply a 2D image
        std::vector<int> naxes;
        naxes.push_back(m_num_x);
        naxes.push_back(m_num_y);
        if (m_num_maps > 1) {
            for (int i = 0; i < ndim(); ++i) {
                naxes.push_back(m_shape[i]);
            }
        }

        // Allocate image
        hdu = new GFitsImageDouble(naxes);

        // Store data in image
        if (naxes.size() == 2) {
            const double* ptr = m_pixels.data();
            for (int iy = 0; iy < m_num_y; ++iy) {
                for (int ix = 0; ix < m_num_x; ++ix) {
                    (*hdu)(ix,iy) = *ptr++;
                }
            }
        }
        else {
            const double* ptr = m_pixels.data();
            for (int imap = 0; imap < m_num_maps; ++imap) {
                for (int iy = 0; iy < m_num_y; ++iy) {
                    for (int ix = 0; ix < m_num_x; ++ix) {
                        (*hdu)(ix,iy,imap) = *ptr++;
                    }
                }
            }
        }

    } // endif: there were pixels

    // ... otherwise create an empty header
    else {
        hdu = new GFitsImageDouble;
    }

    // Set extension name
    hdu->extname("IMAGE");

    // If we have sky projection information then write into FITS header
    if (m_proj != NULL) {
    	m_proj->write(*hdu);
	}

    // Set additional keywords
    //TODO

    // Return HDU
    return hdu;
}


/***********************************************************************//**
 * @brief Compute solid angle subtended by 4 sky directions
 *
 * @param[in] dir1 First sky direction.
 * @param[in] dir2 Second sky direction.
 * @param[in] dir3 Third sky direction.
 * @param[in] dir4 Forth sky direction.
 * @return Solid angle (steradians).
 *
 * Estimate the solid angle subtended by 4 sky directions using Huilier's
 * theorem.
 *
 * Below, the definiton of the pixel cornes and sides are shown as used
 * within the code.
 *
 *             a12
 *         1---------2
 *         |\       /|
 *         | \a13  / |
 *         |  \   /  |
 *         |   \ /   |
 *      a14|    X    |a23
 *         |   / \   |
 *         |  /   \  |
 *         | /a24  \ |
 *         |/       \|
 *         4---------3
 *             a34
 *
 ***************************************************************************/
double GSkyMap::solidangle(const GSkyDir& dir1, const GSkyDir& dir2,
                           const GSkyDir& dir3, const GSkyDir& dir4) const
{
    // Initialise solid angle
    double solidangle = 0.0;

    // Compute angular distances between pixel corners
    double a12 = dir1.dist(dir2);
    double a14 = dir1.dist(dir4);
    double a23 = dir2.dist(dir3);
    double a24 = dir2.dist(dir4);
    double a34 = dir3.dist(dir4);

    // Special case: a12 or a14 is zero, then pixel is a triangle composed
    // of [2,3,4]
    if (a12 <= 0.0 || a14 <= 0.0) {
        double s   = 0.5 * (a23 + a34 + a24);
        solidangle = 4.0 * std::atan(std::sqrt(std::tan(0.5*s) *
                                               std::tan(0.5*(s-a23)) *
                                               std::tan(0.5*(s-a34)) *
                                               std::tan(0.5*(s-a24))));
    }

    // Special case: a23 or a 34 is zero, then pixel is a triangle composed
    // of [1,2,4]
    else if (a23 <= 0.0 || a34 <= 0.0) {
        double s   = 0.5 * (a12 + a24 + a14);
        solidangle = 4.0 * std::atan(std::sqrt(std::tan(0.5*s) *
                                               std::tan(0.5*(s-a12)) *
                                               std::tan(0.5*(s-a24)) *
                                               std::tan(0.5*(s-a14))));
    }

    // Otherwise we have a polygon
    else {

        // Triangle 1 [1,2,4]
        double s1      = 0.5 * (a12 + a24 + a14);
        double excess1 = std::atan(std::sqrt(std::tan(0.5*s1) *
                                             std::tan(0.5*(s1-a12)) *
                                             std::tan(0.5*(s1-a24)) *
                                             std::tan(0.5*(s1-a14))));

        // Triangle 2 [2,3,4]
        double s2      = 0.5 * (a23 + a34 + a24);
        double excess2 = std::atan(std::sqrt(std::tan(0.5*s2) *
                                             std::tan(0.5*(s2-a23)) *
                                             std::tan(0.5*(s2-a34)) *
                                             std::tan(0.5*(s2-a24))));

        // Determine solid angle
        solidangle = 4.0 * (excess1 + excess2);

    } // endif: we had a polynom

    // Return solid angle
    return solidangle;
}


/***********************************************************************//**
 * @brief Compute solid angle subtended by 3 sky directions
 *
 * @param[in] dir1 First sky direction.
 * @param[in] dir2 Second sky direction.
 * @param[in] dir3 Third sky direction.
 * @return Solid angle (steradians).
 *
 * Estimate the solid angle subtended by 3 sky directions using Huilier's
 * theorem.
 *
 * Below, the definiton of the pixel cornes and sides are shown as used
 * within the code.
 *
 *             a12
 *         1---------2
 *         |        / 
 *         |       /  
 *         |      /   
 *         |     /    
 *      a13|    /a23
 *         |   / 
 *         |  / 
 *         | /
 *         |/
 *         3
 *
 ***************************************************************************/
double GSkyMap::solidangle(const GSkyDir& dir1, const GSkyDir& dir2,
                           const GSkyDir& dir3) const
{
    // Compute angular distances between pixel corners
    double a12 = dir1.dist(dir2);
    double a13 = dir1.dist(dir3);
    double a23 = dir2.dist(dir3);

    // Compute solid angle
    double s          = 0.5 * (a12 + a23 + a13);
    double solidangle = 4.0 * std::atan(std::sqrt(std::tan(0.5*s) *
                                                  std::tan(0.5*(s-a12)) *
                                                  std::tan(0.5*(s-a23)) *
                                                  std::tan(0.5*(s-a13))));

    // Return solid angle
    return solidangle;
}


/***********************************************************************//**
 * @brief Checks whether a circular region overlaps with this map
 *
 * @param[in] region Circular region
 * @return True if circular region overlaps with map, false otherwise
 *
 * The check is done by first testing whether the central pointing
 * position of the observation falls into the sky map. If this is false, 
 * then pixel positions that border the sky map are tested for whether or
 * not they fall into the observation region. The positions tested can be
 * visualized as follows, where '*' marks the positions tested.
 *
 *     *   *   *   *   *   *
 *       +---+---+---+---+
 *     * |0,2|1,2|2,2|3,2| *
 *       +---+---+---+---+
 *     * |0,1|1,1|2,1|3,1| *
 *       +---+---+---+---+
 *     * |0,0|1,0|2,0|3,0| *
 *       +---+---+---+---+
 *     *   *   *   *   *   *
 *
 ***************************************************************************/
bool GSkyMap::overlaps_circle(const GSkyRegionCircle& region) const
{
    // Initialise overlap with true
    bool overlap = false;

    // If the centre of the region is inside the map then signal overlap ...
    if (contains(region.centre())) {
        overlap = true;
    }

    // ... otherwise loop through each of the outer bins and check if the
    // region overlaps with them (within some delta)
    else {

        // Initialise test pixels
        GSkyPixel pix1(0.0,0.0);
        GSkyPixel pix2(0.0,0.0);

        // Check the bottom and top bins
        pix1.y(-1.0);
        pix2.y(double(ny()));
        for (double xbin = -1.0; xbin <= nx(); xbin += 1.0) {
            pix1.x(xbin);
            pix2.x(xbin);
            if ((region.contains(pix2dir(pix1))) ||
                (region.contains(pix2dir(pix2)))) {
                overlap = true;
                break;
            }
        } // endfor: loop on xbin

        // If no overlap has been found then check now the left and right
        // bins
        if (!overlap) {
            pix1.x(-1.0);
            pix2.x(double(nx()));
            for (double ybin = 0.0; ybin < ny(); ybin += 1.0) {
                pix1.y(ybin);
                pix2.y(ybin);
                if ((region.contains(pix2dir(pix1))) ||
                    (region.contains(pix2dir(pix2)))) {
                    overlap = true;
                    break;
                }
            } // endfor: loop on ybin
        }

    } // endelse: check border bins

    // Return overlap flag
    return overlap;
}


/***********************************************************************//**
 * @brief Check if HDU contains HEALPix data
 *
 * @param[in] hdu FITS Header Data Unit.
 * @return True is HDU contains HEALPix data.
 *
 * Returns true if the HDU is not an image and the HDU has the "PIXTYPE"
 * keyword set to "HEALPIX".
 ***************************************************************************/
bool GSkyMap::is_healpix(const GFitsHDU& hdu) const
{
    // Initialise flag
    bool flag(false);

    // If PIXTYPE keyword equals "HEALPIX" then signal that we have
    // HEALPix data
    if ((hdu.exttype() != GFitsHDU::HT_IMAGE) &&
        (hdu.has_card("PIXTYPE") && (hdu.string("PIXTYPE") == "HEALPIX"))) {
        flag = true;
    }

    // Return flag
    return (flag);
}


/***********************************************************************//**
 * @brief Check if HDU contains WCS data
 *
 * @param[in] hdu FITS Header Data Unit.
 * @return True is HDU contains WCS data.
 *
 * Returns true if the HDU is an image and the HDU has the "NAXIS" keyword
 * set to a value >= 2.
 ***************************************************************************/
bool GSkyMap::is_wcs(const GFitsHDU& hdu) const
{
    // Initialise flag
    bool flag(false);

    // If extension is an image thn signal that we have WCS data
    if ((hdu.exttype() == GFitsHDU::HT_IMAGE) &&
        (hdu.has_card("NAXIS") && (hdu.integer("NAXIS") >= 2))) {
        flag = true;
    }

    // Return flag
    return (flag);
}


/***********************************************************************//**
 * @brief Check if map is the same
 *
 * @param[in] map Sky map.
 * @return True is sky map definition is identical.
 *
 * Returns true if the sky map definition is identical to the one of the
 * current map.
 ***************************************************************************/
bool GSkyMap::is_same(const GSkyMap& map) const
{
    // Initialise flag
    bool flag = true;

    // Single loop that can be exited any time
    do {

        // If sky projection differs then set flag to false and break
        if (m_proj->code() != map.m_proj->code()) {
            flag = false;
            break;
        }

        // If sky projection has different coordinate system then set flag
        // to false and break
        if (m_proj->coordsys() != map.m_proj->coordsys()) {
            flag = false;
            break;
        }

        // If sky map dimension differs then set flag to false and break
        if ((m_num_x != map.m_num_x) || (m_num_y != map.m_num_y)) {
            flag = false;
            break;
        }

        // For WCS projections, if sky maps have different reference values,
        // pixels or pixel sizes then set flag to false and break
        const GWcs* wcs1 = dynamic_cast<const GWcs*>(m_proj);
        const GWcs* wcs2 = dynamic_cast<const GWcs*>(map.m_proj);
        if ((wcs1 != NULL) && (wcs2 != NULL)) {
            if ((wcs1->crval(0) != wcs2->crval(0)) ||
                (wcs1->crval(1) != wcs2->crval(1)) ||
                (wcs1->crpix(0) != wcs2->crpix(0)) ||
                (wcs1->crpix(1) != wcs2->crpix(1)) ||
                (wcs1->cdelt(0) != wcs2->cdelt(0)) ||
                (wcs1->cdelt(1) != wcs2->cdelt(1))) {
                flag = false;
                break;
            }
        }

        // For Healpix projections, if sky maps have different number of pixels,
        // nside or ordering then set flag to false and break
        const GHealpix* healpix1 = dynamic_cast<const GHealpix*>(m_proj);
        const GHealpix* healpix2 = dynamic_cast<const GHealpix*>(map.m_proj);
        if ((healpix1 != NULL) && (healpix2 != NULL)) {
            if ((healpix1->npix()     != healpix2->npix()) ||
                (healpix1->nside()    != healpix2->nside()) ||
                (healpix1->ordering() != healpix2->ordering())) {
                flag = false;
                break;
            }
        }

    } while(false);

    // Return flag
    return (flag);
}


/***********************************************************************//**
 * @brief Convolve sky map
 *
 * @param[in] kernel Convolution kernel type ("DISK", "GAUSSIAN").
 * @param[in] par Convolution parameter.
 * @param[in] normalise Normalise kernel?
 *
 * Convolves all sky maps using the specified @p kernel and a convolution
 * parameter.
 *
 * For the "DISK" kernel the convolution parameter is the disk radius in
 * degrees. For the "GAUSSIAN" kernel the convolution parameter is the
 * Gaussian sigma in degrees.
 *
 * If @p normalise is true the kernel is normalised to an integral of unity
 * so that the operation corresponds to a smoothing. Otherwise the maximum
 * of the kernel is 1 so that the convolution is a correlation.
 ***************************************************************************/
void GSkyMap::convolve(const std::string& kernel,
                       const double&      par,
                       const bool&        normalise)
{
    // Continue only if the convolution parameter is positive
    if (par > 0.0) {

        // Get FFT of convolution kernel
        GFft fft_kernel(convolution_kernel(kernel, par, normalise));

        // Loop over all sky maps
        for (int k = 0; k < m_num_maps; ++k) {

            // Extract sky map
            GNdarray array(m_num_x, m_num_y);
            const double *src = m_pixels.data() + k*m_num_pixels;
            double       *dst = array.data();
            for (int i = 0; i < m_num_pixels; ++i) {
                *dst++ = *src++;
            }

            // FFT of sky map
            GFft fft_array = GFft(array);

            // Multiply FFT of sky map with FFT of kernel
            GFft fft_smooth = fft_array * fft_kernel;

            // Backward transform sky map
            GNdarray smooth = fft_smooth.backward();

            // Insert sky map
            src = smooth.data();
            dst = m_pixels.data() + k*m_num_pixels;
            for (int i = 0; i < m_num_pixels; ++i) {
                *dst++ = *src++;
            }

        } // endfor: looped over all sky map

    } // endif: convolution parameter was positive

    // Return
    return;
}

/***********************************************************************//**
 * @brief Return convolution kernel
 *
 * @param[in] kernel Convolution kernel type ("DISK", "GAUSSIAN").
 * @param[in] par Convolution parameter (>0).
 * @param[in] normalise Normalise kernel?
 * @return Array filled with convolution kernel.
 *
 * If @p normalise is true the kernel is normalised to an integral of unity
 * so that the operation corresponds to a smoothing. Otherwise the maximum
 * of the kernel is 1 so that the convolution is a correlation.
 ***************************************************************************/
GNdarray GSkyMap::convolution_kernel(const std::string& kernel,
                                     const double&      par,
                                     const bool&        normalise) const
{
    // Get pointer on WCS projection
    const GWcs* wcs = dynamic_cast<const GWcs*>(m_proj);
    if (wcs == NULL) {
        std::string msg = "Sky map is not a WCS projection. Method is only "
                          "valid for WCS projections.";
        throw GException::invalid_argument(G_CONVOLUTION_KERNEL, msg);
    }

    // Get X and Y step size
    double dx = wcs->cdelt(0);
    double dy = wcs->cdelt(1);

    // Initialise kernel
    GNdarray kern(m_num_x, m_num_y);

    // Initialise kernel sum
    double sum = 0.0;

    // Handle disk kernel
    if (gammalib::toupper(kernel) == "DISK") {

        // Fill kernel
        for (int ix1 = 0, ix2 = m_num_x; ix1 < m_num_x; ++ix1, --ix2) {
            double x   = ix1 * dx;
            double xqs = x * x;
            for (int iy1 = 0, iy2 = m_num_y; iy1 < m_num_y; ++iy1, --iy2) {
                double y = iy1 * dy;
                if (std::sqrt(xqs + y*y) <= par) {
                    kern(ix1,iy1) += 1.0;
                    sum           += 1.0;
                    if (ix2 < m_num_x) {
                        kern(ix2,iy1) += 1.0;
                        sum           += 1.0;
                    }
                    if (iy2 < m_num_y) {
                        kern(ix1,iy2) += 1.0;
                        sum           += 1.0;
                    }
                    if ((ix2 < m_num_x) && (iy2 < m_num_y)) {
                        kern(ix2,iy2) += 1.0;
                        sum           += 1.0;
                    }
                }
            }
        }

    } // endif: handled disk kernel

    // Handle Gaussian kernel
    else if (gammalib::toupper(kernel) == "GAUSSIAN") {

        // Compute the exponential normalisation
        double norm = -0.5 / (par * par);

        // Fill kernel
        for (int ix1 = 0, ix2 = m_num_x; ix1 < m_num_x; ++ix1, --ix2) {
            double x   = ix1 * dx;
            double xqs = x * x;
            for (int iy1 = 0, iy2 = m_num_y; iy1 < m_num_y; ++iy1, --iy2) {
                double y       = iy1 * dy;
                double value   = std::exp(norm * (xqs + y*y));
                kern(ix1,iy1) += value;
                sum           += value;
                if (ix2 < m_num_x) {
                    kern(ix2,iy1) += value;
                    sum           += value;
                }
                if (iy2 < m_num_y) {
                    kern(ix1,iy2) += value;
                    sum           += value;
                }
                if ((ix2 < m_num_x) && (iy2 < m_num_y)) {
                    kern(ix2,iy2) += value;
                    sum           += value;
                }
            }
        }

    } // endif: handled Gaussian kernel

    // ... otherwise throw an exception
    else {
        std::string msg = "Invalid kernel type \""+kernel+"\" specified. "
                          "Please specify one of \"DISK\", \"GAUSSIAN\".";
        throw GException::invalid_argument(G_CONVOLUTION_KERNEL, msg);
    }

    // Normalise kernel
    if ((normalise) && (sum > 0.0)) {
        for (int ix = 0; ix < m_num_x; ++ix) {
            for (int iy = 0; iy < m_num_y; ++iy) {
                kern(ix,iy) /= sum;
            }
        }
    }

    // Return kernel
    return kern;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Computes square root of sky map elements
 *
 * @param[in] map Sky map.
 * @return Sky map containing the square root of every element.
 ***************************************************************************/
GSkyMap sqrt(const GSkyMap& map)
{
    // Initialise result vector
    GSkyMap result(map);

    // Loop over all maps
    for (int i = 0; i < map.nmaps(); ++i) {

        // Loop over all bins
        for (int j = 0; j < map.npix(); ++j) {

            // Get the content from the bin
            double content = map(j,i);

            // Throw an exception if content is negative
            if (content < 0.0) {
                std::string msg = "Negative value encountered. "
                                  "Cannot take the sqrt from a negative "
                                  "value";
                throw GException::invalid_value(G_SQRT, msg);
            }

            // Set content of the result map
            result(j,i) = std::sqrt(content);

        } // endfor: looped over all bins

    } // endfor: looped over all maps

    // Return sky map
    return result;
}


/***********************************************************************//**
 * @brief Computes the natural logarithm of sky map elements
 *
 * @param[in] map Sky map.
 * @return Sky map containing the natural logarithm of every element.
 ***************************************************************************/
GSkyMap log(const GSkyMap& map)
{
    // Initialise result vector
    GSkyMap result(map);

    // Loop over all maps
    for (int i = 0; i < map.nmaps(); ++i) {

        // Loop over all bins
        for (int j = 0; j < map.npix(); ++j) {

            // Get the content from the bin
            double content = map(j,i);

            // Throw an exception if content is negative
            if (content <= 0.0) {
                std::string msg = "Negative or zero value encountered. "
                                  "Cannot calculate the logarithm of a "
                                  "negative or zero value";
                throw GException::invalid_value(G_LOG, msg);
            }

            // Set content of the result map
            result(j,i) = std::log(content);

        } // endfor: looped over all bins

    } // endfor: looped over all maps

    // Return sky map
    return result;
}


/***********************************************************************//**
 * @brief Computes the base 10 logarithm of sky map elements
 *
 * @param[in] map Sky map.
 * @return Sky map containing the base 10 logarithm of every element.
 ***************************************************************************/
GSkyMap log10(const GSkyMap& map)
{
    // Initialise result vector
    GSkyMap result(map);

    // Loop over all maps
    for (int i = 0; i < map.nmaps(); ++i) {

        // Loop over all bins
        for (int j = 0; j < map.npix(); ++j) {

            // Get the content from the bin
            double content = map(j,i);

            // Throw an exception if content is negative
            if (content <= 0.0) {
                std::string msg = "Negative or zero value encountered. "
                                  "Cannot calculate the logarithm of a "
                                  "negative or zero value";
                throw GException::invalid_value(G_LOG10, msg);
            }

            // Set content of the result map
            result(j,i) = std::log10(content);

        } // endfor: looped over all bins

    } // endfor: looped over all maps

    // Return sky map
    return result;
}


/***********************************************************************//**
 * @brief Computes the absolute value of sky map elements
 *
 * @param[in] map Sky map.
 * @return Sky map containing the absolute value of every element.
 ***************************************************************************/
GSkyMap abs(const GSkyMap& map)
{
    // Initialise result vector
    GSkyMap result(map);

    // Loop over all maps
    for (int i = 0; i < map.nmaps(); ++i) {

        // Loop over all bins
        for (int j = 0; j < map.npix(); ++j) {

            // Get the content from the bin
            double content = map(j,i);

            // Set content of the result map
            result(j,i) = std::abs(content);

        } // endfor: looped over all bins

    } // endfor: looped over all maps

    // Return sky map
    return result;
}


/***********************************************************************//**
 * @brief Computes the sign value of sky map elements
 *
 * @param[in] map Sky map.
 * @return Sky map containing the sign value of every pixel.
 *
 * This method returns a sky map filled with a value of 1 if the pixel
 * is positive, a value of -1 if the pixel is negative or a value of 0
 * if the pixel is 0.
 ***************************************************************************/
GSkyMap sign(const GSkyMap& map)
{
    // Initialise result vector
    GSkyMap result(map);

    // Loop over all maps
    for (int i = 0; i < map.nmaps(); ++i) {

        // Loop over all bins
        for (int j = 0; j < map.npix(); ++j) {

            // Get the content from the bin
            double content = map(j,i);

            // Handle the 3 cases
            if (content < 0.0) {
                result(j,i) = -1.0;
            }
            else if (content > 0.0) {
                result(j,i) = +1.0;
            }
            else {
                result(j,i) = 0.0;
            }

        } // endfor: looped over all bins

    } // endfor: looed over all maps

    // Return sky map
    return result;
}


/***********************************************************************//**
 * @brief Clips map at given value
 *
 * @param[in] map Sky map.
 * @param[in] thresh Threshold value.
 * @return Sky map containing the value in the original map if >= thresh,
 * otherwise thresh 
 ***************************************************************************/
GSkyMap clip(const GSkyMap& map, const double& thresh)
{
    // Initialise result vector
    GSkyMap result(map);

    // Loop over all maps
    for (int i = 0; i < map.nmaps(); ++i) {

        // Loop over all bins
        for (int j = 0; j < map.npix(); ++j) {

            // Get the content from the bin
            double content = map(j,i);

            // If content is below threshold, set to thresh
            if (content < thresh) {
                result(j,i) = thresh;
            }

            // ... otherwise keep content
            else {
                result(j,i) = content;
            }

        } // endfor: looped over all bins

    } // endfor: looped over all maps

    // Return sky map
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Equality operator
 *
 * @param[in] a First sky map.
 * @param[in] b Second sky map.
 * @return True if @p a and @p b are identical.
 *
 * Two sky maps are considered identical if they have the same projections,
 * coordinate definiton and number of pixels. The actual content of the map
 * does not need to be identical.
 ***************************************************************************/
bool operator==(const GSkyMap &a, const GSkyMap &b)
{
    // Return result
    return a.is_same(b);
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] a First sky map.
 * @param[in] b Second sky map.
 * @return True if @p a and @p b are not identical.
 *
 * Two sky maps are considered different if they either differ in projections,
 * coordinate definiton or number of pixels. The actual content of the map
 * does is not relevant.
 ***************************************************************************/
bool operator!=(const GSkyMap &a, const GSkyMap &b)
{
    // Return result
    return !(a == b);
}
