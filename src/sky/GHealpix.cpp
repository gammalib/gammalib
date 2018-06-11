/***************************************************************************
 *                 GHealpix.cpp - Healpix projection class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2015 by Juergen Knoedlseder                         *
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
 * @file GHealpix.cpp
 * @brief HealPix projection class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GVector.hpp"
#include "GHealpix.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT    "GHealpix::GHealpix(int& ,std::string& ,std::string&)"
#define G_READ                                    "GHealpix::read(GFitsHDU&)"
#define G_XY2DIR                               "GHealpix::xy2dir(GSkyPixel&)"
#define G_DIR2XY2                                "GHealpix::dir2xy(GSkyDir&)"
#define G_NEST2RING                               "GHealpix::nest2ring(int&)"
#define G_RING2NEST                               "GHealpix::ring2nest(int&)"
#define G_PIX2ANG_RING        "GHealpix::pix2ang_ring(int, double*, double*)"
#define G_PIX2ANG_NEST        "GHealpix::pix2ang_nest(int, double*, double*)"
#define G_ORDERING_SET                     "GHealpix::ordering(std::string&)"
#define G_INTERPOLATOR             "GHealpix::interpolator(double&, double&)"
#define G_NEIGHBOURS                       "GHealpix::neighbours(GSkyPixel&)"
#define G_BOUNDARIES                 "GHealpix::boundaries(GSkyPixel&, int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */
const int jrll[12]           = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
const int jpll[12]           = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
const int nb_xoffset[]       = {-1,-1, 0, 1, 1, 1, 0,-1};
const int nb_yoffset[]       = { 0, 1, 1, 1, 0,-1,-1,-1};
const int nb_facearray[][12] = {{  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },  // S
                                {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },  // SE
                                { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },  // E
                                {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },  // SW
                                {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },  // center
                                {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },  // NE
                                { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },  // W
                                {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },  // NW
                                {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 }}; // N
const int nb_swaparray[][3]  = {{ 0,0,3 },   // S
                                { 0,0,6 },   // SE
                                { 0,0,0 },   // E
                                { 0,0,5 },   // SW
                                { 0,0,0 },   // center
                                { 5,0,0 },   // NE
                                { 0,0,0 },   // W
                                { 6,0,0 },   // NW
                                { 3,0,0 }}; // N
const int order_max          = 13;
const int ns_max             = 1 << order_max;

/* __ Static conversion arrays ___________________________________________ */
static short ctab[0x100];
static short utab[0x100];

/* __ Globals ____________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GHealpix::GHealpix(void) : GSkyProjection()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] nside Number of sides.
 * @param[in] order Pixel ordering ('RING' or 'NESTED').
 * @param[in] coords Coordinate system ('EQU' or 'GAL').
 *
 * @exception GException::wcs_hpx_bad_nside 
 *            Invalid nside parameter.
 * @exception GException::wcs_bad_coords 
 *            Invalid coordsys parameter.
 * @exception GException::wcs_hpx_bad_ordering 
 *            Invalid ordering parameter.
 ***************************************************************************/
GHealpix::GHealpix(const int&         nside,
                   const std::string& order,
                   const std::string& coords) : GSkyProjection()
{
    // Initialise class members
    init_members();

    // Check nside parameter (power of 2 between 1 and 8192)
    if (nside2order(nside) == -1) {
        throw GException::wcs_hpx_bad_nside(G_CONSTRUCT, nside);
    }

    // Set coordinate system
    coordsys(coords);

    // Set pixel ordering
    ordering(order);

    // Set Healpix parameters
    m_nside      = nside;
    m_npface     = m_nside * m_nside;
    m_ncap       = 2 * (m_npface - m_nside);
    m_num_pixels = 12 * m_npface;
    m_fact2      = 4.0 / m_num_pixels;
    m_fact1      = 2 * m_nside * m_fact2;
    m_solidangle = gammalib::fourpi / m_num_pixels;
    m_order      = nside2order(m_nside);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor from FITS HDU table
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
GHealpix::GHealpix(const GFitsHDU& hdu) : GSkyProjection()
{
    // Initialise class members
    init_members();

    // Read Healpix definition from FITS table
    read(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] proj Healpix projection.
 ***************************************************************************/
GHealpix::GHealpix(const GHealpix& proj) : GSkyProjection(proj)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(proj);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GHealpix::~GHealpix()
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
 * @param[in] proj Healpix projection.
 * @return Healpix projection.
 ***************************************************************************/
GHealpix& GHealpix::operator=(const GHealpix& proj)
{
    // Execute only if object is not identical
    if (this != &proj) {

        // Copy base class members
        this->GSkyProjection::operator=(proj);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(proj);

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
 * @brief Clear object
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GHealpix::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GSkyProjection::free_members();

    // Initialise members
    this->GSkyProjection::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GHealpix* GHealpix::clone(void) const
{
    return new GHealpix(*this);
}


/***********************************************************************//**
 * @brief Read Healpix definition from FITS header
 *
 * @param[in] hdu FITS HDU.
 *
 * @exception GException::wcs
 *            Unable to load Healpix definition from HDU.
 * @exception GException::wcs_bad_coords
 *            Invalid coordsys parameter.
 * @exception GException::wcs_hpx_bad_ordering
 *            Invalid ordering parameter.
 ***************************************************************************/
void GHealpix::read(const GFitsHDU& hdu)
{
    // Clear object
    clear();

    // Check if we have a healpix representation
    if (hdu.string("PIXTYPE") != "HEALPIX") {
        throw GException::wcs(G_READ, "HDU does not contain Healpix data");
    }

    // Get pixel ordering. First search for the ORDERING keyword, then
    // search for the ORDER keyword
    std::string ordering;
    if (hdu.has_card("ORDERING")) {
        ordering = hdu.string("ORDERING");
    }
    else if (hdu.has_card("ORDER")) {
        ordering = hdu.string("ORDER");
    }

    // Get coordinate system. First search for HIER_CRD keyword (this has been
    // used in older versions of LAT exposure cubes). If not found then search
    // for standard COORDSYS keyword.
    std::string coordsys;
    if (hdu.has_card("HIER_CRD")) {
        coordsys = hdu.string("HIER_CRD");
    }
    else if (hdu.has_card("COORDSYS")) {
        coordsys = hdu.string("COORDSYS");
    }

    // Set coordinate system
    this->coordsys(coordsys);

    // Set pixel ordering
    this->ordering(ordering);

    // Get Healpix resolution and determine number of pixels and solid angle
    m_nside      = hdu.integer("NSIDE");
    m_npface     = m_nside * m_nside;
    m_ncap       = 2 * (m_npface - m_nside);
    m_num_pixels = 12 * m_npface;
    m_fact2      = 4.0 / m_num_pixels;
    m_fact1      = 2 * m_nside * m_fact2;
    m_solidangle = gammalib::fourpi / m_num_pixels;
    m_order      = nside2order(m_nside);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Healpix definition into FITS HDU
 *
 * @param[in] hdu FITS HDU.
 *
 * Writes the following keywords in the FITS HDU:
 * EXTNAME  = HEALPIX
 * PIXTYPE  = HEALPIX
 * NSIDE    = nside() = m_nside
 * FIRSTPIX = 0
 * LASTPIX  = npix()-1 = m_num_pixels-1
 * ORDERING = ordering()
 * COORDSYS = coordsys()
 ***************************************************************************/
void GHealpix::write(GFitsHDU& hdu) const
{
    // Set extension name
    hdu.extname("HEALPIX");

    // Set pixtype (kludge to avoid that a boolean is written)
    std::string pixtype = "HEALPIX";

    // Set keywords
    hdu.card("PIXTYPE",  pixtype,    "HEALPix pixelisation");
    hdu.card("NSIDE",    nside(),    "HEALPix resolution parameter");
    hdu.card("FIRSTPIX", 0,          "Index of first pixel");
    hdu.card("LASTPIX",  npix()-1,   "Index of last pixel");
    hdu.card("ORDERING", ordering(),
             "Pixel ordering scheme, either RING or NESTED");
    hdu.card("COORDSYS", coordsys(),
             "Coordinate system, either EQU or GAL");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pixel Sky map pixel.
 ***************************************************************************/
GSkyDir GHealpix::pix2dir(const GSkyPixel& pixel) const
{
    // Throw an exception if sky map pixel is not 1D
    if (!pixel.is_1D()) {
        std::string msg = "Sky map pixel "+pixel.print()+" is not"
                          " 1-dimensional.\n"
                          "Only 1-dimensional pixels are supported by the"
                          " Healpix projection.";
        throw GException::invalid_argument(G_XY2DIR, msg);
    }

    // Perform ordering dependent conversion
    double theta = 0.0;
    double phi   = 0.0;
    switch (m_ordering) {
    case 0:
        pix2ang_ring(int(pixel), &theta, &phi);
        break;
    case 1:
        pix2ang_nest(int(pixel), &theta, &phi);
        break;
    default:
        break;
    }

    // Store coordinate system-dependent result
    GSkyDir result;
    switch (m_coordsys) {
    case 0:
        result.radec(phi, gammalib::pihalf - theta);
        break;
    case 1:
        result.lb(phi, gammalib::pihalf - theta);
        break;
    default:
        break;
    }
    
    // Return
    return result;
}


/***********************************************************************//**
 * @brief Returns sky map pixel of sky coordinate
 *
 * @param[in] dir Sky coordinate.
 * @return Sky map pixel.
 ***************************************************************************/
GSkyPixel GHealpix::dir2pix(const GSkyDir& dir) const
{
    // Compute coordinate system dependent (z,phi)
    double z = 0;
    double phi = 0;
    switch (m_coordsys) {
    case 0:
        z   = cos(gammalib::pihalf - dir.dec());
        phi = dir.ra();
        break;
    case 1:
        z   = cos(gammalib::pihalf - dir.b());
        phi = dir.l();
        break;
    default:
        break;
    }

    // Perform ordering dependent conversion
    int index;
    switch (m_ordering) {
    case 0:
        index = ang2pix_z_phi_ring(z, phi);
        break;
    case 1:
        index = ang2pix_z_phi_nest(z, phi);
        break;
    default:
        break;
    }

    // Return sky map pixel
    return (GSkyPixel(index));
}


/***********************************************************************//**
 * @brief Return interpolator for given sky direction
 *
 * @param[in] dir Sky direction
 ***************************************************************************/
GBilinear GHealpix::interpolator(const GSkyDir& dir) const
{
    // Compute coordinate system dependent theta and phi (in radians)
    double theta = 0;
    double phi   = 0;
    switch (m_coordsys) {
    case 0:
        theta = gammalib::pihalf - dir.dec();
        phi   = dir.ra();
        break;
    case 1:
        theta = gammalib::pihalf - dir.b();
        phi   = dir.l();
        break;
    default:
        break;
    }

    // Normalize theta and phi value
    theta = gammalib::modulo(theta, gammalib::twopi);
    if (theta > gammalib::pi) {
        phi  += gammalib::pi;
        theta = gammalib::twopi - theta;
    }
    phi = gammalib::modulo(phi, gammalib::twopi);

    // Get interpolator
    GBilinear interpolator = this->interpolator(theta, phi);

    // Return interpolator
    return (interpolator);
}


/***********************************************************************//**
 * @brief Returns ordering parameter.
 ***************************************************************************/
std::string GHealpix::ordering(void) const
{
    // Set pixel ordering type
    std::string ordering;
    switch (m_ordering) {
    case 0:
        ordering = "RING";
        break;
    case 1:
        ordering = "NESTED";
        break;
    default:
        ordering = "UNKNOWN";
        break;
    }

    // Return ordering
    return ordering;
}


/***********************************************************************//**
 * @brief Set pixel ordering.
 *
 * @param[in] ordering Pixel ordering (RING or NEST/NESTED).
 *
 * @exception GException::wcs_hpx_bad_ordering
 *            Invalid ordering parameter.
 ***************************************************************************/
void GHealpix::ordering(const std::string& ordering)
{
    // Convert argument to upper case
    std::string uordering = gammalib::toupper(ordering);

    // Set pixel ordering
    if (uordering == "RING") {
        m_ordering = 0;
    }
    else if (uordering == "NESTED" || uordering == "NEST") {
        m_ordering = 1;
    }
    else {
        throw GException::wcs_hpx_bad_ordering(G_ORDERING_SET, ordering);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return neighbouring pixels of a pixel
 *
 * @param[in] pixel Pixel index.
 * @return Array of neighbouring pixels.
 *
 * Returns the 8 neighbours of a given pixel. The method returns a vector
 * with contains the pixel indices of the SW, W, NW, N, NE, E, SE and S
 * neighbours of @p pix (in the given order). If a neighbour does not exist
 * (this can only be the case for the W, N, E and S neighbors), its entry is
 * set to -1.
 *
 * This method has been adapted from the neighbors() function located in the
 * file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
std::vector<int> GHealpix::neighbours(const GSkyPixel& pixel) const
{
    // Throw an exception if sky map pixel is not 1D
    if (!pixel.is_1D()) {
        std::string msg = "Sky map pixel "+pixel.print()+" is not"
                          " 1-dimensional.\n"
                          "Only 1-dimensional pixels are supported by the"
                          " Healpix projection.";
        throw GException::invalid_argument(G_NEIGHBOURS, msg);
    }

    // Initialise result array
    std::vector<int> result;

    // Determine pixel index and face number
    int ix;
    int iy;
    int face_num;
    pix2xyf(int(pixel), &ix, &iy, &face_num);

    // ...
    const int nsm1 = m_nside - 1;

    // ...
    if ((ix > 0) && (ix < nsm1) && (iy > 0) && (iy < nsm1)) {

        // Ring ordering scheme ...
        if (m_ordering == 0) { // Ring?
            for (int m = 0; m < 8; ++m) {
                int index = xyf2ring(ix+nb_xoffset[m],iy+nb_yoffset[m],face_num);
                result.push_back(index);
            }
        }

        // ... or nested scheme?
        else {
            int fpix = int(face_num) << (2*m_order);
            int px0  = spread_bits(ix);
            int py0  = spread_bits(iy) << 1;
            int pxp  = spread_bits(ix+1);
            int pyp  = spread_bits(iy+1) << 1;
            int pxm  = spread_bits(ix-1);
            int pym  = spread_bits(iy-1) << 1;
            result.push_back(fpix + pxm + py0);
            result.push_back(fpix + pxm + pyp);
            result.push_back(fpix + px0 + pyp);
            result.push_back(fpix + pxp + pyp);
            result.push_back(fpix + pxp + py0);
            result.push_back(fpix + pxp + pym);
            result.push_back(fpix + px0 + pym);
            result.push_back(fpix + pxm + pym);
        }

    } // endif
    
    // ...
    else {
        for (int i = 0; i < 8; ++i) {
            int x     = ix + nb_xoffset[i];
            int y     = iy + nb_yoffset[i];
            int nbnum = 4;
            if (x < 0) {
                x     += m_nside;
                nbnum -= 1;
            }
            else if (x >= m_nside) {
                x     -= m_nside;
                nbnum += 1;
            }
            if (y < 0) {
                y     += m_nside;
                nbnum -= 3;
            }
            else if (y >= m_nside) {
                y     -= m_nside;
                nbnum += 3;
            }

            // Compute face
            int f = nb_facearray[nbnum][face_num];
        
            // If face is valid then compute index
            if (f >= 0) {
                int bits = nb_swaparray[nbnum][face_num>>2];
                if (bits & 1) {
                    x = m_nside - x - 1;
                }
                if (bits & 2) {
                    y = m_nside - y - 1;
                }
                if (bits & 4) {
                    std::swap(x,y);
                }
                if (m_ordering == 0) { // Ring?
                    result.push_back(xyf2ring(x, y, f));
                }
                else {
                    result.push_back(xyf2nest(x, y, f));
                }
            }

            // ... otherwise push back an invalid pixel
            else {
                result.push_back(-1);
            }

        } // endfor
        
    } // endelse

    // Return indices
    return result;
}


/***********************************************************************//**
 * @brief Return pixel boundaries
 *
 * @param[in] pixel Sky pixel.
 * @param[in] step Number of returned points (4*step, defaults to 4).
 * @return Array of pixel boundaries.
 *
 * The method returns a vector with sky directions along the boundary of a
 * HealPix pixel. By default, the 4 corners of HealPix pixel will be
 * returned. The first point corresponds to the northernmost corner, the
 * subsequent points follow the pixel boundary through west, south and east
 * corners. If step>1 more intermediate points will be added after the 4th
 * corner of the pixel.
 *
 * This method has been adapted from the boundaries() function located in the
 * file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
std::vector<GSkyDir> GHealpix::boundaries(const GSkyPixel& pixel,
                                          const int&       step) const
{
    // Throw an exception if sky map pixel is not 1D
    if (!pixel.is_1D()) {
        std::string msg = "Sky map pixel "+pixel.print()+" is not"
                          " 1-dimensional.\n"
                          "Only 1-dimensional pixels are supported by the"
                          " Healpix projection.";
        throw GException::invalid_argument(G_BOUNDARIES, msg);
    }

    // Allocate boundaries
    std::vector<GSkyDir> boundaries;

    // Determine pixel index and face number
    int ix;
    int iy;
    int face;
    pix2xyf(int(pixel), &ix, &iy, &face);

    // ...
    double dc = 0.5        / m_nside;
    double xc = (ix + 0.5) / m_nside;
    double yc = (iy + 0.5) / m_nside;
    double d  = 1.0 / (step*m_nside);
    
    //
    for (int i = 0; i < step; ++i) {
    
        // Declare local coordinates
        double  z;
        double  phi;

        // First coordinate
        xyf2loc(xc+dc-i*d, yc+dc, face, &z, &phi);
        boundaries.push_back(loc2dir(z, phi));

        // Second coordinate
        xyf2loc(xc-dc, yc+dc-i*d, face, &z, &phi);
        boundaries.push_back(loc2dir(z, phi));
                
        // Third coordinate
        xyf2loc(xc-dc+i*d, yc-dc, face, &z, &phi);
        boundaries.push_back(loc2dir(z, phi));
        
        // Forth coordinate
        xyf2loc(xc+dc, yc-dc+i*d, face, &z, &phi);
        boundaries.push_back(loc2dir(z, phi));
        
    } // endfor: looped over step size

    // Return boundaries
    return boundaries;
}


/***********************************************************************//**
 * @brief Return maximum angular distance between pixel centre and corners
 *
 * @return Maximum angular distance between pixel centre and corners (radians).
 *
 * Returns the maximum angular distance (in radians) between any pixel
 * centre and its corners.
 *
 * This method has been adapted from the max_pixrad() function located in the
 * file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
double GHealpix::max_pixrad(void) const
{
    // Compute ...
    double t1 = 1.0 - 1.0/m_nside;
    t1       *= t1;
    
    // Get vectors
    GVector va = set_z_phi(2.0/3.0, gammalib::pi/(4*m_nside));
    GVector vb = set_z_phi(1.0 - t1/3, 0.0);

    // Get angle
    double angle = std::atan2(norm(cross(va, vb)), va * vb);

    // Return angle
    return angle;
}


/***********************************************************************//**
 * @brief Print WCS information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing WCS information.
 ***************************************************************************/
std::string GHealpix::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GHealpix ===");

        // Append information
        result.append("\n"+gammalib::parformat("Coordinate system"));
        result.append(coordsys());
        result.append("\n"+gammalib::parformat("Nside (# of divisions)"));
        result.append(gammalib::str(m_nside));
        result.append("\n"+gammalib::parformat("Npface (pixels per face)"));
        result.append(gammalib::str(m_npface));
        result.append("\n"+gammalib::parformat("Ncap (# of cap pixels)"));
        result.append(gammalib::str(m_ncap));
        result.append("\n"+gammalib::parformat("Npix (# of pixels)"));
        result.append(gammalib::str(m_num_pixels));
        result.append("\n"+gammalib::parformat("Order"));
        result.append(gammalib::str(m_order));
        result.append("\n"+gammalib::parformat("Solid angle per pixel"));
        result.append(gammalib::str(m_solidangle)+" sr");
        result.append("\n"+gammalib::parformat("Ordering")+ordering());

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
void GHealpix::init_members(void)
{
    // Initialise members
    m_nside       = 0;
    m_npface      = 0;
    m_ncap        = 0;
    m_order       = 0;
    m_ordering    = 0;
    m_num_pixels  = 0;
    m_fact1       = 0.0;
    m_fact2       = 0.0;
    m_solidangle  = 0.0;

    // Construct conversion arrays
    for (int m = 0; m < 0x100; ++m) {
    ctab[m] =
         (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
      | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4);
    utab[m] =
         (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
      | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] proj Healpix projection.
 ***************************************************************************/
void GHealpix::copy_members(const GHealpix& proj)
{
    // Copy attributes
    m_nside      = proj.m_nside;
    m_npface     = proj.m_npface;
    m_ncap       = proj.m_ncap;
    m_order      = proj.m_order;
    m_ordering   = proj.m_ordering;
    m_num_pixels = proj.m_num_pixels;
    m_fact1      = proj.m_fact1;
    m_fact2      = proj.m_fact2;
    m_solidangle = proj.m_solidangle;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GHealpix::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns true if argument is identical
 *
 * @param[in] proj Sky projection.
 *
 * This method is a helper for the sky projection friends.
 ***************************************************************************/
bool GHealpix::compare(const GSkyProjection& proj) const
{
    // Initialise result
    bool result = false;
    
    // Continue only we compare to a GHealpix object
    const GHealpix* ptr = dynamic_cast<const GHealpix*>(&proj);
    if (ptr != NULL) {
        result = ((m_coordsys   == ptr->m_coordsys) &&
                  (m_nside      == ptr->m_nside)    &&
                  (m_npface     == ptr->m_npface)   &&
                  (m_ncap       == ptr->m_ncap)     &&
                  (m_order      == ptr->m_order)    &&
                  (m_ordering   == ptr->m_ordering) &&
                  (m_num_pixels == ptr->m_num_pixels));
    }

    // Return result
    return result;
}



/*==========================================================================
 =                                                                         =
 =                         Low-level Healpix methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Compress Bits
 *
 * @param[in] value Value.
 * @return Compressed Bits.
 *
 * This method has been adapted from the compress_bits() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::compress_bits(const int& value) const
{
    // Compress Bits
    int raw        = (value & 0x5555) | ((value & 0x55550000) >> 15);
    int compressed = ctab[raw & 0xff] | (ctab[raw >> 8] << 4);

    // Return compressed value
    return compressed;
}


/***********************************************************************//**
 * @brief Spread Bits
 *
 * @param[in] value Compressed value.
 * @return Spread Bits.
 *
 * This method has been adapted from the spread_bits() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::spread_bits(const int& value) const
{
    // Spread bits
    int spread = utab[value & 0xff] | (utab[(value >> 8) & 0xff] << 16);

    // Return spread value
    return spread;
}


/***********************************************************************//**
 * @brief Convert pixel number in to (x,y,face) tuple
 *
 * @param[in] pix Pixel number.
 * @param[out] ix X index.
 * @param[out] iy Y index.
 * @param[out] face Face number.
 *
 * This method has been adapted from the pix2xyf() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
void GHealpix::pix2xyf(const int& pix, int* ix, int* iy, int* face) const
{
    // Handle ring scheme ...
    if (m_ordering == 0) {
        ring2xyf(pix, ix, iy, face);
    }

    // ... or nested scheme
    else {
        nest2xyf(pix, ix, iy, face);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert pixel number in nested scheme to (x,y,face) tuple
 *
 * @param[in] pix Pixel number in nested scheme.
 * @param[out] ix X index.
 * @param[out] iy Y index.
 * @param[out] face Face number.
 *
 * This method has been adapted from the nest2xyf() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
void GHealpix::nest2xyf(const int& pix, int* ix, int* iy, int* face) const
{
    // Compute face number
    *face = pix >> (2 * m_order);

    // Compute pixel
    int pixel = pix & (m_npface - 1);

    // Compute (x,y)
    *ix = compress_bits(pixel);
    *iy = compress_bits(pixel >> 1);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert pixel number in ring scheme to (x,y,face) tuple
 *
 * @param[in] pix Pixel number in ring scheme
 * @param[out] ix X index.
 * @param[out] iy Y index.
 * @param[out] face Face number.
 *
 * This method has been adapted from the ring2xyf() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
void GHealpix::ring2xyf(const int& pix, int* ix, int* iy, int* face) const
{
    // Declare some variables
    int nl2 = 2*m_nside;
    int iring;
    int iphi;
    int kshift;
    int nr;

    // Handle pixel in the North Polar cap
    if (pix < m_ncap) {
        iring  = (1+isqrt(1+2*pix)) >> 1;    // Counted from North pole
        iphi   = (pix+1) - 2*iring*(iring-1);
        kshift = 0;
        nr     = iring;
        *face  = (iphi-1)/nr;
    }

    // Handle pixel in equatorial region
    else if (pix < (m_num_pixels-m_ncap)) {
        int ip   = pix - m_ncap;
        int tmp  = (m_order>=0) ? ip >> (m_order+2) : ip/(4*m_nside);
        iring    = tmp + m_nside;
        iphi     = ip - tmp * 4 * m_nside + 1;
        kshift   = (iring + m_nside) & 1;
        nr       = m_nside;
        int ire  = iring - m_nside + 1;
        int irm  = nl2 + 2 - ire;
        int ifm  = iphi - ire/2 + m_nside - 1;
        int ifp  = iphi - irm/2 + m_nside - 1;
        if (m_order >= 0) {
            ifm >>= m_order;
            ifp >>= m_order;
        }
        else {
            ifm /= m_nside;
            ifp /= m_nside;
        }
        *face = (ifp==ifm) ? (ifp|4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
    
    // Handle pixel in the South Polar cap
    else {
        int ip = m_num_pixels - pix;
        iring  = (1+isqrt(2*ip-1))>>1; // Counted from South pole
        iphi   = 4 * iring + 1 - (ip - 2*iring*(iring-1));
        kshift = 0;
        nr     = iring;
        iring  = 2 * nl2 - iring;
        *face  = 8 + (iphi-1)/nr;
    }

    // Now compute the (ix,iy) values
    int irt = iring    - (jrll[*face] * m_nside) + 1;
    int ipt = 2 * iphi -  jpll[*face] * nr - kshift -1;
    if (ipt >= nl2) {
        ipt -= 8*m_nside;
    }
    *ix =  (ipt-irt) >> 1;
    *iy = (-ipt-irt) >> 1;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert (x,y,face) tuple to pixel number in nested scheme
 *
 * @param[in] ix X index
 * @param[in] iy Y index
 * @param[in] face Face number
 * @return Pixel number if nested scheme
 *
 * This method has been adapted from the xyf2nest() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::xyf2nest(const int& ix, const int& iy, const int& face) const
{
    // Computed pixel number
    int pix = (int(face) << (2 * m_order)) +
              spread_bits(ix) + (spread_bits(iy) << 1);

    // Return pixel number
    return pix;
}


/***********************************************************************//**
 * @brief Convert (x,y,face) tuple to pixel number in ring scheme
 *
 * @param[in] ix X index
 * @param[in] iy Y index
 * @param[in] face Face number
 * @return Pixel number if ring scheme
 *
 * This method has been adapted from the xyf2ring() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::xyf2ring(const int& ix, const int& iy, const int& face) const
{
    // Compute ring number
    int nl4 = 4 * m_nside;
    int jr  = (jrll[face]*m_nside) - ix - iy  - 1;

    // Get information about that ring
    int  n_before;
    int  nr;
    bool shifted;
    get_ring_info(jr, &n_before, &nr, &shifted);
    
    // Compute pixel number
    nr   >>= 2;
    int kshift = 1-shifted;
    int jp = (jpll[face]*nr + ix - iy + 1 + kshift) / 2;

    //planck_assert(jp<=4*nr,"must not happen");
    
    // Assumption: if this triggers, then nl4==4*nr
    if (jp < 1) {
        jp += nl4;
    }

    // Return pixel number
    return (n_before + jp - 1);
}


/***********************************************************************//**
 * @brief Convert (x,y,f) tuple into local coordinates
 *
 * @param[in] x X value.
 * @param[in] y Y value.
 * @param[in] face Face number.
 * @param[out] z Z value.
 * @param[out] phi Phi value.
 *
 * This method has been adapted from the xyf2loc() function located in the
 * file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
void GHealpix::xyf2loc(const double& x, const double& y, const int& face,
                       double* z, double* phi) const
{
    // ...
    double jr = jrll[face] - x - y;
    double nr;

    // Compute z
    if (jr < 1) {
        nr         = jr;
        double tmp = nr*nr / 3.0;
        *z = 1 - tmp;
    }
    else if (jr > 3) {
        nr         = 4 - jr;
        double tmp = nr*nr / 3.0;
        *z = tmp - 1;
    }
    else {
        nr = 1;
        *z = (2.0-jr) * 2.0/3.0;
    }

    // Compute Phi
    double tmp = jpll[face] * nr + x - y;
    if (tmp < 0) {
        tmp += 8;
    }
    if (tmp >= 8) {
        tmp -= 8;
    }
    *phi = (nr < 1.0e-15) ? 0.0 : (0.5 * gammalib::pihalf * tmp) / nr;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert local coordinate into sky direction
 *
 * @param[in] z Z value.
 * @param[in] phi Phi value.
 * @return Sky direction.
 *
 * This method has been adapted from the locToVec3() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
GSkyDir GHealpix::loc2dir(const double& z, const double& phi) const
{
    // Compute longitude and latitude
    double sintheta  = std::sqrt((1.0 - z) * (1.0 + z));
    double longitude = std::atan2(sintheta * std::sin(phi), sintheta * std::cos(phi));
    double latitude  = std::asin(z);

    // Set sky direction
    GSkyDir dir;
    switch (m_coordsys) {
    case 0:
        dir.radec(longitude, latitude);
        break;
    case 1:
        dir.lb(longitude, latitude);
        break;
    default:
        break;
    }

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Convert nside to order
 *
 * @param[in] nside Number of sides.
 * @return Order of HealPix projection.
 ***************************************************************************/
int GHealpix::nside2order(const int& nside) const
{
    // Initialise order
    int order = -1;

    // Determine order
    for (int m = 0; m <= order_max; ++m) {
        int nstest = 1 << m;
        if (nside == nstest) {
            order = m;
            break;
        }
        if (nside < nstest)
            break;
    }

    // Return order
    return order;
}


/***********************************************************************//**
 * @brief Converts pixel number in nested indexing scheme to ring scheme
 *
 * @param[in] pix Pixel number in nested indexing scheme.
 * @return Pixel number in ring indexing scheme.
 *
 * @exception GException::invalid_value
 *            Healpix projection does not represent a hiearachical map.
 *
 * This method has been adapted from the nest2ring() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::nest2ring(const int& pix) const
{
    // Throw an exception if map is not a hierachical map
    if (m_order < 0) {
        std::string msg = "A hierarchical map projection is required.";
        throw GException::invalid_value(G_NEST2RING, msg);
    }

    // Convert nested index to (x,y,face) tuple
    int ix;
    int iy;
    int face;
    nest2xyf(pix, &ix, &iy, &face);

    // Convert (x,y,face) tuple to ring index
    int iring = xyf2ring(ix, iy, face);

    // Return ring index
    return iring;
}


/***********************************************************************//**
 * @brief Converts pixel number in ring indexing scheme to nested scheme
 *
 * @param[in] pix Pixel number in ring indexing scheme.
 * @return Pixel number in nested indexing scheme.
 *
 * @exception GException::invalid_value
 *            Healpix projection does not represent a hiearachical map.
 *
 * This method has been adapted from the ring2nest() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::ring2nest(const int& pix) const
{
    // Throw an exception if map is not a hierachical map
    if (m_order < 0) {
        std::string msg = "A hierarchical map projection is required.";
        throw GException::invalid_value(G_RING2NEST, msg);
    }

    // Convert nested index to (x,y,face) tuple
    int ix;
    int iy;
    int face;
    ring2xyf(pix, &ix, &iy, &face);

    // Convert (x,y,face) tuple to ring index
    int iring = xyf2nest(ix, iy, face);

    // Return ring index
    return iring;
}


/***********************************************************************//**
 * @brief Convert pixel index to (x,y) coordinate
 *
 * @param[in] ipix Pixel index for which (x,y) are to be computed.
 * @param[out] x Pointer to x coordinate.
 * @param[out] y Pointer to y coordinate.
 ***************************************************************************/
void GHealpix::pix2xy(const int& ipix, int* x, int* y) const
{
    // Set x coordinate
    int raw = (ipix & 0x5555) | ((ipix & 0x55550000) >> 15);
    *x      = ctab[raw & 0xff] | (ctab[raw >> 8] << 4);

    // Set y coordinate
    raw = ((ipix & 0xaaaa) >> 1) | ((ipix & 0xaaaa0000) >> 16);
    *y  = ctab[raw & 0xff] | (ctab[raw >> 8] << 4);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert (x,y) coordinate to pixel
 *
 * @param[in] x x coordinate.
 * @param[in] y y coordinate.
 ***************************************************************************/
int GHealpix::xy2pix(int x, int y) const
{
    // Return pixel
    return utab[x&0xff] | (utab[x>>8]<<16) | (utab[y&0xff]<<1) | (utab[y>>8]<<17);
}


/***********************************************************************//**
 * @brief Convert pixel index to (theta,phi) angles for ring ordering
 *
 * @param[in] ipix Pixel index for which (theta,phi) are to be computed.
 * @param[out] theta Pointer to result zenith angle in radians.
 * @param[out] phi Pointer to result azimuth angle in radians.
 *
 * @exception GException::out_of_range
 *            Pixel index is out of range.
 ***************************************************************************/
void GHealpix::pix2ang_ring(int ipix, double* theta, double* phi) const
{
    // Check if ipix is in range
    if (ipix < 0 || ipix >= m_num_pixels) {
        throw  GException::out_of_range(G_PIX2ANG_RING, ipix, 0, m_num_pixels-1);
    }

    // Handle North Polar cap
    if (ipix < m_ncap) {
        int iring = int(0.5*(1+isqrt(1+2*ipix))); // counted from North pole
        int iphi  = (ipix+1) - 2*iring*(iring-1);
        *theta    = std::acos(1.0 - (iring*iring) * m_fact2);
        *phi      = (iphi - 0.5) * gammalib::pi/(2.0*iring);
    }

    // Handle Equatorial region
    else if (ipix < (m_num_pixels - m_ncap)) {
        int    ip    = ipix - m_ncap;
        int    iring = ip/(4*m_nside) + m_nside;   // counted from North pole
        int    iphi  = ip%(4*m_nside) + 1;
        double fodd  = ((iring+m_nside)&1) ? 1 : 0.5;
        int    nl2   = 2*m_nside;
        *theta       = std::acos((nl2 - iring) * m_fact1);
        *phi         = (iphi - fodd) * gammalib::pi/nl2;
    }

    // Handle South Polar cap
    else {
        int ip    = m_num_pixels - ipix;
        int iring = int(0.5*(1+isqrt(2*ip-1)));    // Counted from South pole
        int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
        *theta    = std::acos(-1.0 + (iring*iring) * m_fact2);
        *phi      = (iphi - 0.5) * gammalib::pi/(2.*iring);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert pixel index to (theta,phi) angles for nested ordering
 *
 * @param[in] ipix Pixel index for which (theta,phi) are to be computed.
 * @param[out] theta Pointer to result zenith angle in radians.
 * @param[out] phi Pointer to result azimuth angle in radians.
 *
 * @exception GException::out_of_range
 *            Pixel index is out of range.
 ***************************************************************************/
void GHealpix::pix2ang_nest(int ipix, double* theta, double* phi) const
{
    // Check if ipix is in range
    if (ipix < 0 || ipix >= m_num_pixels) {
        throw GException::out_of_range(G_PIX2ANG_NEST, ipix, 0, m_num_pixels-1);
    }

    // Get face number and index in face
    int nl4      = 4 * m_nside;
    int face_num = ipix >> (2*m_order);      // Face number in {0,11}
    int ipf      = ipix & (m_npface - 1);

    // Get pixel coordinates
    int ix;
    int iy;
    pix2xy(ipf, &ix, &iy);

    // Computes the z coordinate on the sphere
    int jr = (jrll[face_num] << m_order) - ix - iy - 1;

    // Declare result variables
    int    nr;
    double z;
    int    kshift;

    // North pole region
    if (jr < m_nside) {
        nr     = jr;
        z      = 1. - nr*nr*m_fact2;
        kshift = 0;
    }

    // South pole region
    else if (jr > 3*m_nside) {
        nr     = nl4 - jr;
        z      = nr*nr*m_fact2 - 1;
        kshift = 0;
    }

    // Equatorial region
    else {
        nr     = m_nside;
        z      = (2*m_nside-jr) * m_fact1;
        kshift = (jr-m_nside) & 1;
    }

    // Computes the phi coordinate on the sphere, in [0,2Pi]
    int jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
    if (jp > nl4) jp -= nl4;
    if (jp <   1) jp += nl4;

    // Computes Theta and Phi
    *theta = std::acos(z);
    *phi   = (jp - (kshift+1)*0.5) * (gammalib::pihalf / nr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pixel which contains angular coordinates (z,phi)
 *
 * @param[in] z Cosine of zenith angle - cos(theta).
 * @param[in] phi Azimuth angle in radians.
 ***************************************************************************/
int GHealpix::ang2pix_z_phi_ring(double z, double phi) const
{
    // Initialise pixel
    int ipix = 0;

    // Setup
    double za = fabs(z);
    double tt = gammalib::modulo(phi, gammalib::twopi) * gammalib::inv_pihalf; // in [0,4)

    // Equatorial region
    if (za <= gammalib::twothird) {
        double temp1  = m_nside*(0.5+tt);
        double temp2  = m_nside*z*0.75;
        int    jp     = int(temp1-temp2);           // index of ascending edge line
        int    jm     = int(temp1+temp2);           // index of descending edge line
        int    ir     = m_nside + 1 + jp - jm;      // in {1,2n+1}
        int    kshift = 1 - (ir & 1);               // kshift=1 if ir even, 0 otherwise
        int    ip     = (jp+jm-m_nside+kshift+1)/2; // in {0,4n-1}
        ip            = int(gammalib::modulo(ip,4*m_nside));
        ipix          = m_ncap + (ir-1)*4*m_nside + ip;
    }

    // North & South polar caps
    else {
        double tp  = tt - int(tt);
        double tmp = m_nside * std::sqrt(3*(1-za));
        int    jp  = int(tp*tmp);       // increasing edge line index
        int    jm  = int((1.0-tp)*tmp); // decreasing edge line index
        int    ir  = jp + jm + 1;       // ring number counted from the closest pole
        int    ip  = int(tt*ir);        // in {0,4*ir-1}
        ip = int(gammalib::modulo(ip,4*ir));
        if (z>0)
            ipix = 2*ir*(ir-1) + ip;
        else
            ipix = m_num_pixels - 2*ir*(ir+1) + ip;
    }

    // Return pixel
    return ipix;
}


/***********************************************************************//**
 * @brief Returns pixel which contains angular coordinates (z,phi)
 *
 * @param[in] z Cosine of zenith angle - cos(theta).
 * @param[in] phi Azimuth angle in radians.
 ***************************************************************************/
int GHealpix::ang2pix_z_phi_nest(double z, double phi) const
{
    // Initialise face and pixel numbers
    int face_num;
    int ix;
    int iy;

    // Setup
    double za = fabs(z);
    double tt = gammalib::modulo(phi, gammalib::twopi) * gammalib::inv_pihalf; // in [0,4)

    // Equatorial region
    if (za <= gammalib::twothird) {
        double temp1 = ns_max*(0.5+tt);
        double temp2 = ns_max*z*0.75;
        int    jp    = int(temp1-temp2); // index of  ascending edge line
        int    jm    = int(temp1+temp2); // index of descending edge line
        int    ifp   = jp >> order_max;  // in {0,4}
        int    ifm   = jm >> order_max;
        if (ifp == ifm)                  // faces 4 to 7
            face_num = (ifp==4) ? 4: ifp+4;
        else if (ifp < ifm)              // (half-)faces 0 to 3
            face_num = ifp;
        else                             // (half-)faces 8 to 11
            face_num = ifm + 8;
        ix = jm & (ns_max-1);
        iy = ns_max - (jp & (ns_max-1)) - 1;
    }

    // Polar region, za > 2/3
    else {
        int    ntt = int(tt);
        double tp  = tt-ntt;
        double tmp = ns_max * std::sqrt(3*(1-za));
        int    jp  = int(tp*tmp);        // increasing edge line index
        int    jm  = int((1.0-tp)*tmp);  // decreasing edge line index
        if (jp >= ns_max) jp = ns_max-1; // for points too close to the boundary
        if (jm >= ns_max) jm = ns_max-1;
        if (z >= 0) {
            face_num = ntt;              // in {0,3}
            ix       = ns_max - jm - 1;
            iy       = ns_max - jp - 1;
        }
        else {
            face_num = ntt + 8;          // in {8,11}
            ix       =  jp;
            iy       =  jm;
        }
    }

    // Get pixel
    int ipf = xy2pix(ix, iy);
    ipf >>= (2*(order_max - m_order));     // in {0, nside**2 - 1}
    return ipf + (face_num<<(2*m_order));  // in {0, 12*nside**2 - 1}
}


/***********************************************************************//**
 * @brief Get ring index north of cos(theta)
 *
 * @param[in] costheta Cosine of colatitude theta
 *
 * Returns the number of the next ring to the north of @a cos(theta). It may
 * return 0; in this case @a costheta lies north of all rings.
 *
 * This method has been adapted from the ring_above() function located in
 * the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
int GHealpix::ring_above(const double& costheta) const
{
    // Initialise ring index
    int iring;

    // Get absolute cosine
    double acostheta = std::abs(costheta);

    // Are we in the equatorial region
    if (acostheta <= gammalib::twothird) {
        iring = int(m_nside * (2.0 - 1.5 * costheta));
    }

    // ... otherwise we're in the pole region
    else {
        iring = int(m_nside * std::sqrt(3.0 * (1.0 - acostheta)));
        if (costheta <= 0.0) {
            iring = 4 * m_nside - iring - 1;
        }
    }

    // Return
    return iring;
}


/***********************************************************************//**
 * @brief Returns useful information about a given ring of the projection
 *
 * @param[in] ring The ring number (the number of the first ring is 1)
 * @param[out] startpix The number of the first pixel in the ring
 * @param[out] ringpix The number of pixels in the ring
 * @param[out] shifted If @a true, the center of the first pixel is not at
 *                     @a phi=0
 *
 * This method has been adapted from the get_ring_info_small() function
 * located in the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
void GHealpix::get_ring_info(const int& ring,
                             int*       startpix,
                             int*       ringpix,
                             bool*      shifted) const
{
    // Handle ring in North polar cap
    if (ring < m_nside) {
        *shifted  = true;
        *ringpix  = 4 * ring;
        *startpix = 2 * ring * (ring-1);
    }
    
    // Handle ring in equatorial region
    else if (ring < 3*m_nside) {
        *shifted  = ((ring-m_nside) & 1) == 0;
        *ringpix  = 4 * m_nside;
        *startpix = m_ncap + (ring-m_nside) * *ringpix;
    }

    // Handle ring in South polar cap
    else {
        int nr    = 4 * m_nside - ring;
        *shifted  = true;
        *ringpix  = 4 * nr;
        *startpix = m_num_pixels - 2 * nr * (nr+1);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns useful information about a given ring of the projection
 *
 * @param[in] ring The ring number (the number of the first ring is 1)
 * @param[out] startpix The number of the first pixel in the ring
 * @param[out] ringpix The number of pixels in the ring
 * @param[out] theta (radians)
 * @param[out] shifted If @a true, the center of the first pixel is not at
 *                     @a phi=0
 *
 * This method has been adapted from the get_ring_info2() function
 * located in the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
void GHealpix::get_ring_info(const int& ring,
                             int*       startpix,
                             int*       ringpix,
                             double*    theta,
                             bool*      shifted) const
{
    //
    int northring = (ring > 2 * m_nside) ? 4 * m_nside - ring : ring;

    // Are we in the North?
    if (northring < m_nside) {
        double tmp      = northring * northring * m_fact2;
        double costheta = 1.0 - tmp;
        double sintheta = std::sqrt(tmp * (2.0-tmp));
        *startpix       = 2 * northring * (northring - 1);
        *ringpix        = 4 * northring;
        *theta          = std::atan2(sintheta, costheta);
        *shifted        = true;
    }
    else {
        *theta    = std::acos((2.0 * m_nside-northring) * m_fact1);
        *ringpix  = 4 * m_nside;
        *shifted  = ((northring - m_nside) & 1) == 0;
        *startpix = m_ncap + (northring - m_nside) * *ringpix;
    }
    
    // Are we in the southern hemisphere?
    if (northring != ring) {
        *theta    = gammalib::pi - *theta;
        *startpix = m_num_pixels - *startpix - *ringpix;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return interpolator
 *
 * @param[in] theta Colatitude of direction (radian, the North pole is at theta=0)
 * @param[in] phi Longitude of direction (radians)
 *
 * @exception GException::invalid_argument
 *            Invalid @p theta argument
 *
 * Returns a bilinear pixel interpolator for a given sky direction.
 *
 * This method has been adapted from the get_interpol() function
 * located in the file healpix_base.cc in Healpix version 3.20.
 ***************************************************************************/
GBilinear GHealpix::interpolator(const double& theta, const double& phi) const
{
    // Allocate interpolator
    GBilinear interpolator;

    // Check that theta is valid
    if (theta < 0.0 || theta > gammalib::pi) {
        std::string msg = "Colatitude "+gammalib::str(theta)+" is outside "
                          "valid range [0,pi].";
        throw GException::invalid_argument(G_INTERPOLATOR, msg);
    }

    // Prepare computation
    double costheta = std::cos(theta);
    int    ir1      = ring_above(costheta); // Ring above actual colatitude
    int    ir2      = ir1 + 1;              // Ring below actual colatitude
    int    sp;                              // Start pixel in ring
    int    nr;                              // Number of pixels in ring
    double theta1;                          // Colatitude of ring above
    double theta2;                          // Colatitude of ring below
    bool   shift;

    // Compute interpolating pixels and phi weights if the colatitude is not
    // North of all rings (if we're North of all rings, ir1=0)
    if (ir1 > 0) {
        get_ring_info(ir1, &sp, &nr, &theta1, &shift);
        double dphi = gammalib::twopi / nr;  // Phi spacing of pixels
        double tmp  = (phi/dphi - 0.5*shift);
        int    i1   = (tmp < 0.0) ? int(tmp)-1 : int(tmp);
        double w1   = (phi - (i1+0.5*shift)*dphi) / dphi;
        int    i2   = i1 + 1;
        if (i1 < 0) {
            i1 += nr;
        }
        if (i2 >= nr) {
            i2 -= nr;
        }
        interpolator.index1()  = sp + i1;
        interpolator.index2()  = sp + i2;
        interpolator.weight1() = 1.0 - w1;
        interpolator.weight2() = w1;
    }
    
    // Compute interpolating pixels and phi weights is the colatitude is not
    // South of all rings (if we're South of all rings, ir2=4*m_nside)
    if (ir2 < (4*m_nside)) {
        get_ring_info(ir2, &sp, &nr, &theta2, &shift);
        double dphi = gammalib::twopi / nr;  // Phi spacing of pixels
        double tmp  = (phi/dphi - 0.5*shift);
        int    i1   = (tmp < 0.0) ? int(tmp)-1 : int(tmp);
        double w1   = (phi - (i1+0.5*shift)*dphi) / dphi;
        int    i2   = i1 + 1;
        if (i1 < 0) {
            i1 += nr;
        }
        if (i2 >= nr) {
            i2 -= nr;
        }
        interpolator.index3()  = sp + i1;
        interpolator.index4()  = sp + i2;
        interpolator.weight3() = 1.0 - w1;
        interpolator.weight4() = w1;
    }

    // Now handle the special case that the colatitude is  North of all
    // rings
    if (ir1 == 0) {
        double wtheta           = theta/theta2;
        interpolator.weight3() *= wtheta;
        interpolator.weight4() *= wtheta;
        double fac              = (1.0-wtheta)*0.25;
        interpolator.weight1()  = fac;
        interpolator.weight2()  = fac;
        interpolator.weight3() += fac;
        interpolator.weight4() += fac;
        interpolator.index1()   = (interpolator.index3() + 2) & 3;
        interpolator.index2()   = (interpolator.index4() + 2) & 3;
    }

    // ... and now the case that the colatitude is South of all rings
    else if (ir2 == 4*m_nside) {
        double wtheta            = (theta-theta1) / (gammalib::pi-theta1);
        interpolator.weight1()  *= (1.0 - wtheta);
        interpolator.weight2()  *= (1.0 - wtheta);
        double fac              = wtheta*0.25;
        interpolator.weight1() += fac;
        interpolator.weight2() += fac;
        interpolator.weight3()  = fac;
        interpolator.weight4()  = fac;
        interpolator.index3()   = ((interpolator.index1() + 2) & 3) + m_num_pixels - 4;
        interpolator.index4()   = ((interpolator.index2() + 2) & 3) + m_num_pixels - 4;
    }
    
    // ... and now multiply-in the theta weights for the general case
    else {
        double wtheta           = (theta-theta1) / (theta2-theta1);
        interpolator.weight1() *= (1.0 - wtheta);
        interpolator.weight2() *= (1.0 - wtheta);
        interpolator.weight3() *= wtheta;
        interpolator.weight4() *= wtheta;
    }

    // If we have a nested pixel scheme then convert now the ring
    // indices into nested indices
    if (m_ordering == 1) {
        interpolator.index1() = ring2nest(interpolator.index1());
        interpolator.index2() = ring2nest(interpolator.index2());
        interpolator.index3() = ring2nest(interpolator.index3());
        interpolator.index4() = ring2nest(interpolator.index4());
    }

    // Return interpolator
    return interpolator;
}


/***********************************************************************//**
 * @brief Integer n that fulfills n*n <= arg < (n+1)*(n+1)
 *
 * @param[in] arg Argument.
 *
 * Returns the integer @a n, which fulfills @a n*n <= arg < (n+1)*(n+1).
 ***************************************************************************/
unsigned int GHealpix::isqrt(unsigned int arg) const
{
    // Return
    return unsigned(std::sqrt(arg+0.5));
}


/***********************************************************************//**
 * @brief Return 3D vector
 *
 * @param[in] z Z value.
 * @param[in] phi Phi value.
 * @return 3D vector
 *
 * This method has been adapted from the set_z_phi() function located in the
 * file vec3.h in Healpix version 3.20.
 ***************************************************************************/
GVector GHealpix::set_z_phi(const double& z, const double& phi) const
{
    // Initialise 3D vector
    GVector vector(3);

    // Assign elements
    double sintheta = std::sqrt((1.0 - z) * (1.0 + z));
    vector[0] = sintheta * std::cos(phi);
    vector[1] = sintheta * std::sin(phi);
    vector[2] = z;

    // Return vector
    return vector;
}
