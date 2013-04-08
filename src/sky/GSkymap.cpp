/***************************************************************************
 *              GSkymap.cpp - Class that implements a sky map              *
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
 * @file GSkymap.cpp
 * @brief Sky map class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GSkymap.hpp"
#include "GWcsRegistry.hpp"
#include "GWcslib.hpp"
#include "GWcsHPX.hpp"
#include "GFits.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsImageDouble.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT_HPX     "GSkymap::GSkymap(std::string,std::string,int," \
                                                           "std::string,int)"
#define G_CONSTRUCT_MAP "GSkymap::GSkymap(std::string,std::string,GSkyDir," \
                                                 "int,int,double,double,int)"
#define G_OP_ACCESS_1D                         "GSkymap::operator(int&,int&)"
#define G_OP_ACCESS_2D                   "GSkymap::operator(GSkyPixel&,int&)"
#define G_OP_VALUE                         "GSkymap::operator(GSkyDir&,int&)"
#define G_READ                               "GSkymap::read(const GFitsHDU*)"
#define G_PIX2DIR                                     "GSkymap::pix2dir(int)"
#define G_DIR2PIX                                 "GSkymap::dir2pix(GSkyDir)"
#define G_XY2DIR                                 "GSkymap::xy2dir(GSkyPixel)"
#define G_DIR2XY                                   "GSkymap::dir2xy(GSkyDir)"
#define G_OMEGA1                                        "GSkymap::omega(int)"
#define G_OMEGA2                                  "GSkymap::omega(GSkyPixel)"
#define G_SET_WCS "GSkymap::set_wcs(std::string,std::string,double,double," \
                               "double,double,double,double,GMatrix,GVector)"
#define G_READ_HEALPIX                   "GSkymap::read_healpix(GFitsTable*)"
#define G_READ_WCS                           "GSkymap::read_wcs(GFitsImage*)"
#define G_ALLOC_WCS                         "GSkymap::alloc_wcs(GFitsImage*)"

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
 ***************************************************************************/
GSkymap::GSkymap(void)
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
 ***************************************************************************/
GSkymap::GSkymap(const std::string& filename)
{
    // Initialise class members for clean destruction
    init_members();

    // Load skymap
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Healpix constructor
 *
 * @param[in] wcs World Coordinate System (HPX).
 * @param[in] coords Coordinate System (CEL or GAL).
 * @param[in] nside Nside parameter.
 * @param[in] order Pixel ordering (RING or NEST).
 * @param[in] nmaps Number of maps in set (default=1).
 *
 * @exception GException::wcs_invalid
 *            Invalid wcs parameter.
 * @exception GException::skymap_bad_par
 *            Invalid sky map parameter.
 ***************************************************************************/
GSkymap::GSkymap(const std::string& wcs, const std::string& coords,
                 const int& nside, const std::string& order,
                 const int nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    // Check if wcs is HPX
    if (toupper(wcs) != "HPX") {
        throw GException::wcs_invalid(G_CONSTRUCT_HPX, wcs,
                                      "WCS parameter must be 'HPX'.");
    }

    // Check if nmaps parameter is >0
    if (nmaps < 1) {
        throw GException::skymap_bad_par(G_CONSTRUCT_HPX, nmaps,
                                         "nmaps parameter must be >0.");
    }

    // Allocate WCS
    m_wcs = new GWcsHPX(nside, order, coords);

    // Set number of pixels and number of maps
    m_num_pixels = static_cast<GWcsHPX*>(m_wcs)->npix();
    m_num_maps   = nmaps;

    // Allocate pixels
    alloc_pixels();

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
 * @exception GException::skymap_bad_par
 *            Invalid sky map parameter.
 ***************************************************************************/
GSkymap::GSkymap(const std::string& wcs, const std::string& coords,
                 double const& x, double const& y,
                 double const& dx, double const& dy,
                 const int& nx, const int& ny, const int nmaps)
{
    // Initialise class members for clean destruction
    init_members();

    // Check parameters
    if (nx < 1) {
        throw GException::skymap_bad_par(G_CONSTRUCT_MAP, nx,
                                         "nx parameter must be >0.");
    }
    if (ny < 1) {
        throw GException::skymap_bad_par(G_CONSTRUCT_MAP, ny,
                                         "ny parameter must be >0.");
    }
    if (nmaps < 1) {
        throw GException::skymap_bad_par(G_CONSTRUCT_MAP, nmaps,
                                         "nmaps parameter must be >0.");
    }

    // Set WCS
    double  crval1 = x;
    double  crval2 = y;
    double  crpix1 = double(nx+1)/2.0;
    double  crpix2 = double(ny+1)/2.0;
    double  cdelt1 = dx;
    double  cdelt2 = dy;
    GMatrix cd(2,2);
    GVector pv2(21);
    set_wcs(wcs, coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2, 
            cd, pv2);

    // Set number of pixels and number of maps
    m_num_x      = nx;
    m_num_y      = ny;
    m_num_pixels = m_num_x * m_num_y;
    m_num_maps   = nmaps;

    // Allocate pixels
    alloc_pixels();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] map Sky map.
 ***************************************************************************/
GSkymap::GSkymap(const GSkymap& map)
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
GSkymap::~GSkymap(void)
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
 * @brief 1D pixel access operator
 *
 * @param[in] pixel Pixel index (0,1,...,m_num_pixels).
 * @param[in] map Map index (0,1,...,m_num_maps).
 *
 * @exception GException::out_of_range
 *            Pixel index and/or map index are outside valid range.
 *
 * Access sky map pixel by its index, where the most quickly varying axis is
 * the x axis of the map.
 *
 * @todo Implement proper skymap exception (actual is for matrix elements)
 ***************************************************************************/
double& GSkymap::operator() (const int& pixel, const int& map)
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (pixel < 0 || pixel >= m_num_pixels) {
        throw GException::out_of_range(G_OP_ACCESS_1D, pixel, 0, m_num_pixels-1);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_1D, map, 0, m_num_maps-1);
    }
    #endif

    // Return reference to pixel value
    return m_pixels[pixel+m_num_pixels*map];
}


/***********************************************************************//**
 * @brief 1D pixel access operator (const variant)
 *
 * @param[in] pixel Pixel index (0,1,...,m_num_pixels).
 * @param[in] map Map index (0,1,...,m_num_maps).
 *
 * @exception GException::out_of_range
 *            Pixel index and/or map index are outside valid range.
 *
 * Access sky map pixel by its index, where the most quickly varying axis is
 * the x axis of the map.
 *
 * @todo Implement proper skymap exception (actual is for matrix elements)
 ***************************************************************************/
const double& GSkymap::operator() (const int& pixel, const int& map) const
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (pixel < 0 || pixel >= m_num_pixels) {
        throw GException::out_of_range(G_OP_ACCESS_1D, pixel, 0, m_num_pixels-1);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_1D, map, 0, m_num_maps-1);
    }
    #endif

    // Return reference to pixel value
    return m_pixels[pixel+m_num_pixels*map];
}


/***********************************************************************//**
 * @brief 2D pixel access operator
 *
 * @param[in] pixel Sky pixel.
 * @param[in] map Map index (0,1,...,m_num_maps).
 *
 * @exception GException::out_of_range
 *            Sky pixel and/or map index are outside valid range.
 *
 * Access sky map pixel by its 2D index (x,y) that is implemented by the
 * GSkyPixel class.
 *
 * @todo Implement proper skymap exception (actual is for matrix elements)
 ***************************************************************************/
double& GSkymap::operator() (const GSkyPixel& pixel, const int& map)
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (!isinmap(pixel)) {
        throw GException::out_of_range(G_OP_ACCESS_2D,
                                       int(pixel.x()), int(pixel.y()),
                                       m_num_x-1, m_num_y-1);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_2D, map, 0, m_num_maps-1);
    }
    #endif

    // Get pixel index
    int index = xy2pix(pixel);

    // Return reference to pixel value
    return m_pixels[index+m_num_pixels*map];
}


/***********************************************************************//**
 * @brief 2D pixel access operator
 *
 * @param[in] pixel Sky pixel.
 * @param[in] map Map index (0,1,...,m_num_maps).
 *
 * @exception GException::out_of_range
 *            Sky pixel and/or map index are outside valid range.
 *
 * Access sky map pixel by its 2D index (x,y) that is implemented by the
 * GSkyPixel class.
 *
 * @todo Implement proper skymap exception (actual is for matrix elements)
 ***************************************************************************/
const double& GSkymap::operator() (const GSkyPixel& pixel, const int& map) const
{
    // Throw an error if pixel index or map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (!isinmap(pixel)) {
        throw GException::out_of_range(G_OP_ACCESS_2D,
                                       int(pixel.x()), int(pixel.y()),
                                       m_num_x-1, m_num_y-1);
    }
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_ACCESS_2D, map, 0, m_num_maps-1);
    }
    #endif

    // Get pixel index
    int index = xy2pix(pixel);

    // Return reference to pixel value
    return m_pixels[index+m_num_pixels*map];
}


/***********************************************************************//**
 * @brief Return interpolated skymap value for sky direction
 *
 * @param[in] dir Sky direction.
 * @param[in] map Map index (0,1,...,m_num_maps).
 *
 * @exception GException::out_of_range
 *            Map index lies outside valid range.
 *
 * Returns the skymap value for a given sky direction, obtained by bi-linear
 * interpolation of the neighbouring pixels. If the sky direction falls
 * outside the area covered by the skymap, a value of 0 is returned.
 ***************************************************************************/
double GSkymap::operator() (const GSkyDir& dir, const int& map) const
{
    // Throw an error if the map index is not in valid range
    #if defined(G_RANGE_CHECK)
    if (map < 0 || map >= m_num_maps) {
        throw GException::out_of_range(G_OP_VALUE, map, 0, m_num_maps-1);
    }
    #endif

    // Initialise intensity
    double intensity = 0.0;

    // Determine sky pixel
    GSkyPixel pixel = dir2xy(dir);

    // Continue only if pixel is within the map
    if (isinmap(pixel)) {

        // Set left indices for interpolation. The left index is comprised
        // between 0 and npixels-2. By definition, the right index is then
        // the left index + 1
        int inx_x = int(pixel.x());
        int inx_y = int(pixel.y());
        if (inx_x < 0) {
            inx_x = 0;
        }
        else if (inx_x > m_num_x-2) {
            inx_x = m_num_x - 2;
        }
        if (inx_y < 0) {
            inx_y = 0;
        }
        else if (inx_y > m_num_y-2) {
            inx_y = m_num_y - 2;
        }

        // Set weighting factors for interpolation
        double wgt_x_right = (pixel.x() - inx_x);
        double wgt_x_left  = 1.0 - wgt_x_right;
        double wgt_y_right = (pixel.y() - inx_y);
        double wgt_y_left  = 1.0 - wgt_y_right;

        // Compute skymap pixel indices for bi-linear interpolation
        int inx1 = inx_x + inx_y * m_num_x;
        int inx2 = inx1 + m_num_x;
        int inx3 = inx1 + 1;
        int inx4 = inx2 + 1;

        // Compute weighting factors for bi-linear interpolation
        double wgt1 = wgt_x_left  * wgt_y_left;
        double wgt2 = wgt_x_left  * wgt_y_right;
        double wgt3 = wgt_x_right * wgt_y_left;
        double wgt4 = wgt_x_right * wgt_y_right;

        // Compute map offset
        int offset = m_num_pixels * map;

        // Compute interpolated skymap value
        intensity = wgt1 * m_pixels[inx1 + offset] +
                    wgt2 * m_pixels[inx2 + offset] +
                    wgt3 * m_pixels[inx3 + offset] +
                    wgt4 * m_pixels[inx4 + offset];

    } // endif: pixel was within map

    // Return intensity
    return intensity;
}


/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] map Sky map.
 ***************************************************************************/
GSkymap& GSkymap::operator= (const GSkymap& map)
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


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

 /***********************************************************************//**
 * @brief Clear instance.
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GSkymap::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GSkymap* GSkymap::clone(void) const
{
    return new GSkymap(*this);
}



/***********************************************************************//**
 * @brief Load skymap from FITS file.
 *
 * @param[in] filename FITS file name..
 *
 * Loads HEALPix and non HEALPix skymaps. First searches for HEALPix map in
 * FITS file by scanning all HDUs for PIXTYPE=HEALPIX. If no HEALPix map has
 * been found then search load first non-empty image.
 ***************************************************************************/
void GSkymap::load(const std::string& filename)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Open FITS file
    GFits fits(filename);

    // Get number of HDUs
    int num = fits.size();

    // Initialize pointer to HDU and load flag
    GFitsHDU* hdu    = NULL;
    bool      loaded = false;

    // First search for HEALPix extension. We can skip the first extension
    // since this is always an image and a HEALPix map is stored in a
    // binary table
    for (int extno = 1; extno < num; ++extno) {

        // Get pointer to HDU
        hdu = fits.hdu(extno);

        // If PIXTYPE keyword equals "HEALPIX" then load map
        try {
            if (hdu->string("PIXTYPE") == "HEALPIX") {
                read_healpix(static_cast<const GFitsBinTable*>(hdu));
                loaded = true;
                break;
            }
        }
        catch (GException::fits_key_not_found &e) {
        }

    } // endfor: looped over HDUs

    // If we have not found a HEALPIX map then search now for image.
    // Skip empty images
    if (!loaded) {
        for (int extno = 0; extno < num; ++extno) {

            // Get pointer to HDU
            hdu = fits.hdu(extno);

            // Skip if extension is not an image
            if (extno > 0) {
                if (hdu->string("XTENSION") != "IMAGE")
                    continue;
            }

            // Load WCS map
            read_wcs(static_cast<const GFitsImageDouble*>(hdu));
            loaded = true;
            break;

        } // endfor: looped over HDUs
    } // endif: no HEALPix map found

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save skymap into FITS file.
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * The method does nothing if the skymap holds no valid WCS.
 ***************************************************************************/
void GSkymap::save(const std::string& filename, bool clobber) const
{
    // Continue only if we have data to save
    if (m_wcs != NULL) {

        // Initialise HDU pointer
        GFitsHDU* hdu = NULL;

        // Case A: Skymap is Healpix
        if (m_wcs->code() == "HPX") {
            hdu = create_healpix_hdu();
        }

        // Case B: Skymap is not Healpix
        else {
            hdu = create_wcs_hdu();
        }

        // Create FITS file and save it to disk
        if (hdu != NULL) {
            GFits fits;
            fits.append(*hdu);
            fits.saveto(filename, clobber);
        }

        // Delete HDU
        if (hdu != NULL) delete hdu;

    } // endif: we had data to save

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read skymap from FITS HDU
 *
 * @param[in] hdu FITS HDU.
 *
 * The method returns an empty skymap of the HDU pointer was not valid.
 ***************************************************************************/
void GSkymap::read(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Continue only if HDU pointer is valid
    if (hdu != NULL) {

        // Initialize load flag
        bool loaded = false;

        // Try load as HEALPix map
        try {
            if (hdu->string("PIXTYPE") == "HEALPIX") {
                read_healpix(static_cast<const GFitsBinTable*>(hdu));
                loaded = true;
            }
        }
        catch (GException::fits_key_not_found &e) {
        }

        // ... otherwise try loading as non HEALPix map
        if (!loaded) {

            // Load only if HDU contains an image
            if (hdu->exttype() == 0) {
                read_wcs(static_cast<const GFitsImageDouble*>(hdu));
                loaded = true;
            }

        } // endif

    } // endif: HDU pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write skymap into FITS file
 *
 * @param[in] file FITS file pointer.
 ***************************************************************************/
void GSkymap::write(GFits* file) const
{
    // Continue only if we have data to save
    if (m_wcs != NULL) {

        // Initialise HDU pointer
        GFitsHDU* hdu = NULL;

        // Case A: Skymap is Healpix
        if (m_wcs->code() == "HPX") {
            hdu = create_healpix_hdu();
        }

        // Case B: Skymap is not Healpix
        else {
            hdu = create_wcs_hdu();
        }

        // Append HDU to FITS file.
        if (hdu != NULL) {
            file->append(*hdu);
        }

        // Delete HDU
        if (hdu != NULL) delete hdu;

    } // endif: we had data to save

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 *
 * @exception GException::wcs
 *            No valid WCS found.
 *
 * Returns sky direction for a given sky map pixel. This methods works for
 * sky maps that have a 1D (e.g. HEALPix) or 2D pixel indexation scheme.
 * Automatic index conversion is provided so that 1D schemes need only to
 * implement 1D conversion methods while 2D schemes need only to implement
 * 2D conversion methods.
 ***************************************************************************/
GSkyDir GSkymap::pix2dir(const int& pix) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL) {
        throw GException::wcs(G_PIX2DIR, "No valid WCS found.");
    }

    // Determine sky direction from pixel. Use 2D version if sky map is
    // 2D, otherwise use 1D version.
    GSkyDir dir = (m_num_x == 0) ? m_wcs->pix2dir(pix)
                                 : m_wcs->xy2dir(pix2xy(pix));

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Returns pixel index for a given sky direction
 *
 * @param[in] dir Sky direction.
 *
 * @exception GException::wcs
 *            No valid WCS found.
 *
 * Returns sky map pixel for a given sky direction. This methods works for
 * sky maps that have a 1D (e.g. HEALPix) or 2D pixel indexation scheme.
 * Automatic index conversion is provided so that 1D schemes need only to
 * implement 1D conversion methods while 2D schemes need only to implement
 * 2D conversion methods.
 ***************************************************************************/
int GSkymap::dir2pix(const GSkyDir& dir) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL) {
        throw GException::wcs(G_DIR2PIX, "No valid WCS found.");
    }

    // Determine 1D pixel index for a given sky direction
    int pix = (m_num_x == 0) ? m_wcs->dir2pix(dir)
                             : xy2pix(m_wcs->dir2xy(dir));

    // Return pixel index
    return pix;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Sky map pixel.
 *
 * @exception GException::wcs
 *            No valid WCS found.
 *
 * Returns sky direction for a given sky map pixel. This methods works for
 * sky maps that have a 1D (e.g. HEALPix) or 2D pixel indexation scheme.
 * Automatic index conversion is provided so that 1D schemes need only to
 * implement 1D conversion methods while 2D schemes need only to implement
 * 2D conversion methods.
 ***************************************************************************/
GSkyDir GSkymap::xy2dir(const GSkyPixel& pix) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL) {
        throw GException::wcs(G_XY2DIR, "No valid WCS found.");
    }

    // Determine sky direction from pixel. Use 2D version if sky map is
    // 2D, otherwise use 1D version.
    GSkyDir dir = (m_num_x == 0) ? m_wcs->pix2dir(xy2pix(pix))
                                 : m_wcs->xy2dir(pix);

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Returns sky map pixel for a given sky direction
 *
 * @param[in] dir Sky direction.
 *
 * @exception GException::wcs
 *            No valid WCS found.
 *
 * Returns sky map pixel for a given sky direction. This methods works for
 * sky maps that have a 1D (e.g. HEALPix) or 2D pixel indexation scheme.
 * Automatic index conversion is provided so that 1D schemes need only to
 * implement 1D conversion methods while 2D schemes need only to implement
 * 2D conversion methods.
 ***************************************************************************/
GSkyPixel GSkymap::dir2xy(const GSkyDir& dir) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL) {
        throw GException::wcs(G_DIR2XY, "No valid WCS found.");
    }

    // Determine 1D pixel index for a given sky direction
    GSkyPixel pix = (m_num_x == 0) ? pix2xy(m_wcs->dir2pix(dir))
                                   : m_wcs->dir2xy(dir);

    // Return pixel index
    return pix;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 *
 * @exception GException::wcs
 *            No valid WCS found.
 ***************************************************************************/
double GSkymap::omega(const int& pix) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL) {
        throw GException::wcs(G_OMEGA1, "No valid WCS found.");
    }

    // Determine solid angle from pixel. Use 2D version if sky map is
    // 2D, otherwise use 1D version.
    double omega = (m_num_x == 0) ? m_wcs->omega(pix)
                                  : m_wcs->omega(pix2xy(pix));

    // Return solid angle
    return omega;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Sky map pixel.
 *
 * @exception GException::wcs
 *            No valid WCS found.
 ***************************************************************************/
double GSkymap::omega(const GSkyPixel& pix) const
{
    // Throw error if WCS is not valid
    if (m_wcs == NULL) {
        throw GException::wcs(G_OMEGA2, "No valid WCS found.");
    }

    // Determine solid angle from pixel. Use 2D version if sky map is
    // 2D, otherwise use 1D version.
    double omega = (m_num_x == 0) ? m_wcs->omega(xy2pix(pix))
                                  : m_wcs->omega(pix);

    // Return solid angle
    return omega;
}


/***********************************************************************//**
 * @brief Returns number of pixels
 ***************************************************************************/
int GSkymap::npix(void) const
{
    // Return number of pixels
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns number of pixels in x coordinate
 ***************************************************************************/
int GSkymap::nx(void) const
{
    // Return number of pixels
    return m_num_x;
}


/***********************************************************************//**
 * @brief Returns number of pixels in y coordinate
 ***************************************************************************/
int GSkymap::ny(void) const
{
    // Return number of pixels
    return m_num_y;
}


/***********************************************************************//**
 * @brief Returns number of maps
 ***************************************************************************/
int GSkymap::nmaps(void) const
{
    // Return number of maps
    return m_num_maps;
}


/***********************************************************************//**
 * @brief Converts 2D index (x,y) into 1D pixel index
 *
 * @param[in] pix 2D pixel index.
 *
 * The (x,y) value is rounded to nearest integers before conversion. The x 
 * axis is assumed as the most rapidely varying index.
 ***************************************************************************/
int GSkymap::xy2pix(const GSkyPixel& pix) const
{
    // Get x and y indices by rounding the (x,y) values
    int ix = int(pix.x()+0.5);
    int iy = int(pix.y()+0.5);

    // Return index
    return (ix+iy*m_num_x);
}


/***********************************************************************//**
 * @brief Converts 1D pixel index into 2D index (x,y)
 *
 * @param[in] pix 1D pixel index.
 *
 * If skymap is not 2D the pixel index is returned in the x element of the
 * GSkyPixel object.
 ***************************************************************************/
GSkyPixel GSkymap::pix2xy(const int& pix) const
{
    // Get x and y indices
    double x;
    double y;
    if (m_num_x != 0) {
        x = double(pix % m_num_x);
        y = double(pix / m_num_x);
    }
    else {
        x = double(pix);
        y = 0.0;
    }

    // Set pixel
    GSkyPixel pixel(x, y);

    // Return pixel
    return pixel;
}


/***********************************************************************//**
 * @brief Set WCS skymap
 *
 * @param[in] wcs World Coordinate System.
 *
 * Sets the World Coordinate System of sky map. The method performs a deep
 * copy of the input argument, allowing to destroy the argument after using
 * the method. 
 ***************************************************************************/
void GSkymap::wcs(const GWcs& wcs)
{
    // Free any existing WCS
    if (m_wcs != NULL) delete m_wcs;

    // Clone input WCS
    m_wcs = wcs.clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Verifies if sky direction falls in map
 *
 * @param[in] dir Sky direction.
 *
 * This method checks if the specified sky direction falls within the pixels
 * covered by the skymap. The method uses the dir2xy method to convert the
 * sky direction into 2D pixel indices, and then checks whether the pixel
 * indices fall in the skymap.
 ***************************************************************************/
bool GSkymap::isinmap(const GSkyDir& dir) const
{
    // Convert sky direction into sky pixel
    GSkyPixel pixel = dir2xy(dir);
    
    // Return location flag
    return (isinmap(pixel));
}


/***********************************************************************//**
 * @brief Verifies if sky pixel falls in map
 *
 * @param[in] pixel Sky pixel.
 *
 * This method checks if the specified sky pixel is within the skymap.
 ***************************************************************************/
bool GSkymap::isinmap(const GSkyPixel& pixel) const
{
    // Initialise location flag
    bool inmap = false;

    // If pixel is in range then set location flag to true
    if ((pixel.x()+0.5 >= 0.0 && pixel.x()-0.5 < m_num_x) &&
        (pixel.y()+0.5 >= 0.0 && pixel.y()-0.5 < m_num_y)) {
        inmap = true;
    }

    // Return location flag
    return inmap;
}


/***********************************************************************//**
 * @brief Print sky map
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing sky map information.
 ***************************************************************************/
std::string GSkymap::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSkymap ===");

        // Append information
        result.append("\n"+parformat("Number of pixels")+str(m_num_pixels));
        result.append("\n"+parformat("Number of maps")+str(m_num_maps));
        if (m_wcs != NULL && m_wcs->code() != "HPX") {
            result.append("\n"+parformat("X axis dimension")+str(m_num_x));
            result.append("\n"+parformat("Y axis dimension")+str(m_num_y));
        }

        // Append WCS information
        if (m_wcs != NULL) {
            result.append("\n"+m_wcs->print(chatter));
        }
        else {
            result.append("\n"+parformat("WCS")+"not defined");
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
void GSkymap::init_members(void)
{
    // Initialise members
    m_num_pixels = 0;
    m_num_maps   = 0;
    m_num_x      = 0;
    m_num_y      = 0;
    m_wcs        = NULL;
    m_pixels     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate skymap pixels
 ***************************************************************************/
void GSkymap::alloc_pixels(void)
{
    // Compute data size
    int size = m_num_pixels * m_num_maps;

    // Continue only if there are pixels
    if (size > 0) {

        // Allocate pixels and initialize them to 0
        m_pixels = new double[size];
        for (int i = 0; i < size; ++i) {
            m_pixels[i] = 0.0;
        }

    } // endif: there were pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] map Sky map.
 ***************************************************************************/
void GSkymap::copy_members(const GSkymap& map)
{
    // Copy attributes
    m_num_pixels = map.m_num_pixels;
    m_num_maps   = map.m_num_maps;
    m_num_x      = map.m_num_x;
    m_num_y      = map.m_num_y;

    // Clone WCS if it is valid
    if (map.m_wcs != NULL) m_wcs = map.m_wcs->clone();

    // Compute data size
    int size = m_num_pixels * m_num_maps;

    // Copy pixels
    if (size > 0 && map.m_pixels != NULL) {
        alloc_pixels();
        for (int i = 0; i <  size; ++i) {
            m_pixels[i] = map.m_pixels[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkymap::free_members(void)
{
    // Free memory
    if (m_wcs    != NULL) delete m_wcs;
    if (m_pixels != NULL) delete [] m_pixels;

    // Signal free pointers
    m_wcs        = NULL;
    m_pixels     = NULL;

    // Reset number of pixels
    m_num_pixels = 0;
    m_num_maps   = 0;
    m_num_x      = 0;
    m_num_y      = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set WCS
 *
 * @param[in] wcs World Coordinate System code.
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel.
 * @param[in] crpix2 Y index of reference pixel.
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 * @param[in] cd Astrometry parameters (2x2 matrix, deg/pixel).
 * @param[in] pv2 Projection parameters (length WCS type dependent).
 *
 * @exception GException::wcs_invalid
 *            Invalid wcs parameter (World Coordinate System not known).
 *
 * This method sets the WCS projection pointer based on the WCS code and
 * sky map parameters. It makes use of the GWcsRegistry class to allocate
 * the correct derived class. Note that this method does not support the
 * HPX projection.
 *
 * @todo Remove cd and pv2 parameters.
 ***************************************************************************/
void GSkymap::set_wcs(const std::string& wcs, const std::string& coords,
                      const double& crval1, const double& crval2,
                      const double& crpix1, const double& crpix2,
                      const double& cdelt1, const double& cdelt2,
                      const GMatrix& cd, const GVector& pv2)
{
    // Convert WCS to upper case
    std::string uwcs = toupper(wcs);
    
    // Check if HPX was requested (since this is not allowed)
    if (uwcs == "HPX") {
       throw GException::wcs_invalid(G_SET_WCS, uwcs,
                                     "Method not valid for HPX projection.");
    }
    
    // ... otherwise get projection from registry
    else {
        // Allocate WCS registry
        GWcsRegistry registry;
        
        // Allocate projection from registry
        m_wcs = registry.alloc(uwcs);
        
        // Signal if projection type is not known
        if (m_wcs == NULL) {
            std::string message = "Projection code not known. "
                                  "Should be one of "+registry.list()+".";
            throw GException::wcs_invalid(G_SET_WCS, uwcs, message);
        }

        // Setup WCS
        static_cast<GWcslib*>(m_wcs)->set(coords, crval1, crval2, crpix1, crpix2,
                                          cdelt1, cdelt2);
    
    } // endelse: got projection from registry

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read Healpix data from FITS table.
 *
 * @param[in] hdu FITS HDU containing the Healpix data.
 *
 * HEALPix data may be stored in various formats depending on the 
 * application that has writted the data. HEALPix IDL, for example, may
 * store the data in vectors of length 1024 if the number of pixels is
 * a multiple of 1024. On the other hand, vectors may also be used to store
 * several HEALPix maps into a single column. Alternatively, multiple maps
 * may be stored in multiple columns.
 ***************************************************************************/
void GSkymap::read_healpix(const GFitsTable* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Determine number of rows and columns in table
        int nrows = hdu->nrows();
        int ncols = hdu->ncols();
        #if defined(G_READ_HEALPIX_DEBUG)
        std::cout << "nrows=" << nrows << " ncols=" << ncols << std::endl;
        #endif

        // Allocate Healpix WCS
        m_wcs = new GWcsHPX;

        // Read WCS information from FITS header
        m_wcs->read(hdu);

        // Set number of pixels based on NSIDE parameter
        m_num_pixels = static_cast<GWcsHPX*>(m_wcs)->npix();
        #if defined(G_READ_HEALPIX_DEBUG)
        std::cout << "m_num_pixels=" << m_num_pixels << std::endl;
        #endif

        // Number of map pixels has to be a multiple of the number of
        // rows in column
        if (m_num_pixels % nrows != 0) {
            throw GException::skymap_bad_size(G_READ_HEALPIX, nrows,
                                              m_num_pixels);
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
            const GFitsTableCol* col = &(*hdu)[icol];
            if (col->number() % nentry == 0) {
                m_num_maps += col->number() / nentry;
            }
        }
        #if defined(G_READ_HEALPIX_DEBUG)
        std::cout << "m_num_maps=" << m_num_maps << std::endl;
        #endif

        // Allocate pixels to hold the map
        alloc_pixels();

        // Initialise map counter
        int imap = 0;

        // Loop over all columns
        for (int icol = 0; icol < ncols; ++icol) {

            // Get next column
            const GFitsTableCol* col = &(*hdu)[icol];

            // Only consider columns that can fully hold maps
            if (col->number() % nentry == 0) {

                // Determine number of maps in column
                int num = col->number() / nentry;

                // Loop over all maps in column
                int inx_start = 0;
                int inx_end   = nentry;
                for (int i = 0; i < num; ++i) {

                    // Load map
                    double *ptr = m_pixels + m_num_pixels*imap;
                    for (int row = 0; row < col->length(); ++row) {
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

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read WCS image from FITS HDU
 *
 * @param[in] hdu FITS HDU containing the WCS image.
 *
 * @exception GException::skymap_bad_image_dim
 *            WCS image has invalid dimension (naxis=2 or 3).
 ***************************************************************************/
void GSkymap::read_wcs(const GFitsImage* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Allocate WCS
        alloc_wcs(hdu);

        // Read WCS information from FITS header
        m_wcs->read(hdu);

        // Extract map dimension and number of maps from image
        if (hdu->naxis() == 2) {
            m_num_x    = hdu->naxes(0);
            m_num_y    = hdu->naxes(1);
            m_num_maps = 1;
        }
        else if (hdu->naxis() >= 3) {
            m_num_x    = hdu->naxes(0);
            m_num_y    = hdu->naxes(1);
            m_num_maps = hdu->naxes(2);
        }
        else {
            throw GException::skymap_bad_image_dim(G_READ_WCS, hdu->naxis());
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
        alloc_pixels();

        // Read image
        if (hdu->naxis() == 2) {
            double* ptr = m_pixels;
            for (int iy = 0; iy < m_num_y; ++iy) {
                for (int ix = 0; ix < m_num_x; ++ix) {
                    *ptr++ = hdu->pixel(ix,iy);
                }
            }
        }
        else {
            double* ptr = m_pixels;
            for (int imap = 0; imap < m_num_maps; ++imap) {
                for (int iy = 0; iy < m_num_y; ++iy) {
                    for (int ix = 0; ix < m_num_x; ++ix) {
                        *ptr++ = hdu->pixel(ix,iy,imap);
                    }
                }
            }
        }

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate WCS class
 *
 * @param[in] hdu FITS HDU containing the WCS image.
 *
 * @exception GException::fits_key_not_found
 *            Unable to find required FITS header keyword.
 * @exception GException::skymap_bad_ctype
 *            CTYPE1 and CTYPE2 keywords are incompatible.
 * @exception GException::wcs_invalid
 *            WCS projection of FITS file not supported by GammaLib.
 ***************************************************************************/
void GSkymap::alloc_wcs(const GFitsImage* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get standard keywords
        std::string ctype1 = hdu->string("CTYPE1");
        std::string ctype2 = hdu->string("CTYPE2");

        // Extract projection type
        std::string xproj = ctype1.substr(5,3);
        std::string yproj = ctype2.substr(5,3);

        // Check that projection type is identical on both axes
        if (xproj != yproj) {
            throw GException::skymap_bad_ctype(G_ALLOC_WCS,
                                               ctype1, ctype2);
        }

        // Allocate WCS registry
        GWcsRegistry registry;
        
        // Allocate projection from registry
        m_wcs = registry.alloc(xproj);
        
        // Signal if projection type is not known
        if (m_wcs == NULL) {
            std::string message = "Projection code not known. "
                                  "Should be one of "+registry.list()+".";
            throw GException::wcs_invalid(G_ALLOC_WCS, xproj, message);
        }

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Create FITS HDU containing Healpix data
 *
 * This method allocates a binary table HDU that contains the Healpix data.
 * Deallocation of the table has to be done by the client.
 ***************************************************************************/
GFitsBinTable* GSkymap::create_healpix_hdu(void) const
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
        double* ptr = m_pixels;
        for (int inx = 0; inx < number; ++inx) {
            for (int row = 0; row < rows; ++row) {
                column(row,inx) = *ptr++;
            }
        }

        // Create HDU that contains Healpix map in a binary table
        hdu = new GFitsBinTable(rows);
        hdu->append_column(column);

    } // endif: there were pixels

    // ... otherwise create an empty header
    else {
        hdu = new GFitsBinTable;
    }

    // Set extension name
    hdu->extname("HEALPIX");

    // If we have WCS information then write into FITS header
    if (m_wcs != NULL) m_wcs->write(hdu);

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
GFitsImageDouble* GSkymap::create_wcs_hdu(void) const
{
    // Initialise result to NULL pointer
    GFitsImageDouble* hdu = NULL;

    // Compute size of Healpix data
    int size = m_num_pixels * m_num_maps;

    // Continue only if we have pixels
    if (size > 0) {

        // Set axis parameters for image construction
        int naxis   = (m_num_maps == 1) ? 2 : 3;
        int naxes[] = {m_num_x, m_num_y, m_num_maps};

        // Allocate image
        hdu = new GFitsImageDouble(naxis, naxes);

        // Store data in image
        if (naxis == 2) {
            double* ptr = m_pixels;
            for (int iy = 0; iy < m_num_y; ++iy) {
                for (int ix = 0; ix < m_num_x; ++ix) {
                    (*hdu)(ix,iy) = *ptr++;
                }
            }
        }
        else {
            double* ptr = m_pixels;
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

    // If we have WCS information then write into FITS header
    if (m_wcs != NULL) m_wcs->write(hdu);

    // Set additional keywords
    //TODO

    // Return HDU
    return hdu;
}
