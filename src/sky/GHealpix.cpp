/***************************************************************************
 *                 GHealpix.cpp - Healpix projection class                 *
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
#include "GHealpix.hpp"
//#include "GWcsRegistry.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT           "GHealpix::GHealpix(int,std::string,std::string)"
#define G_READ                                     "GHealpix::read(GFitsHDU*)"
#define G_XY2DIR                                "GHealpix::xy2dir(GSkyPixel&)"
#define G_DIR2XY2                                 "GHealpix::dir2xy(GSkyDir&)"
#define G_PIX2ANG_RING           "GHealpix::pix2ang_ring(int,double*,double*)"
#define G_PIX2ANG_NEST           "GHealpix::pix2ang_nest(int,double*,double*)"
#define G_ORDERING_SET                       "GHealpix::coordsys(std::string)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */
const int jrll[12]  = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
const int jpll[12]  = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
const int order_max = 13;
const int ns_max    = 1 << order_max;

/* __ Static conversion arrays ___________________________________________ */
static short ctab[0x100];
static short utab[0x100];

/* __ Globals ____________________________________________________________ */
//const GHealpix      g_wcs_hpx_seed;
//const GWcsRegistry g_wcs_hpx_registry(&g_wcs_hpx_seed);


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
    m_omega      = gammalib::fourpi / m_num_pixels;
    m_order      = nside2order(m_nside);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor from FITS HDU table
 *
 * @param[in] hdu Pointer to FITS HDU.
 ***************************************************************************/
GHealpix::GHealpix(const GFitsHDU* hdu) : GSkyProjection()
{
    // Initialise class members
    init_members();

    // Read Healpix defintion from FITS HDU
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
 * @brief Read Healpix definiton from FITS header
 *
 * @param[in] hdu FITS HDU containing the Healpix definition.
 *
 * @exception GException::wcs
 *            Unable to load Healpix definition from HDU.
 * @exception GException::wcs_bad_coords
 *            Invalid coordsys parameter.
 * @exception GException::wcs_hpx_bad_ordering
 *            Invalid ordering parameter.
 ***************************************************************************/
void GHealpix::read(const GFitsHDU* hdu)
{
    // Clear object
    clear();

    // Check if we have a healpix representation
    if (hdu->string("PIXTYPE") != "HEALPIX")
        throw GException::wcs(G_READ, "HDU does not contain Healpix data");

    // Get pixel ordering
    std::string order = hdu->string("ORDERING");

    // Get coordinate system.
    // First search for HIER_CRD keyword (this has been used in older
    // versions of LAT exposure cubes). If not found then search for standard
    // COORDSYS keyword.
    std::string coords;
    try {
        coords = hdu->string("HIER_CRD");
    }
    catch (GException::fits_key_not_found &e) {
        coords = hdu->string("COORDSYS");
    }

    // Set coordinate system
    coordsys(coords);

    // Set pixel ordering
    ordering(order);

    // Get Healpix resolution and determine number of pixels and solid angle
    m_nside      = hdu->integer("NSIDE");
    m_npface     = m_nside * m_nside;
    m_ncap       = 2 * (m_npface - m_nside);
    m_num_pixels = 12 * m_npface;
    m_fact2      = 4.0 / m_num_pixels;
    m_fact1      = 2 * m_nside * m_fact2;
    m_omega      = gammalib::fourpi / m_num_pixels;
    m_order      = nside2order(m_nside);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Healpix definiton into FITS HDU
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
void GHealpix::write(GFitsHDU* hdu) const
{
    // Continue only if file pointer is valid
    if (hdu != NULL) {

        // Set extension name
        hdu->extname("HEALPIX");

        // Set keywords
        hdu->card("PIXTYPE",  "HEALPIX",  "HEALPix pixelisation");
        hdu->card("NSIDE",    nside(),    "HEALPix resolution parameter");
        hdu->card("FIRSTPIX", 0,          "Index of first pixel");
        hdu->card("LASTPIX",  npix()-1,   "Index of last pixel");
        hdu->card("ORDERING", ordering(),
                  "Pixel ordering scheme, either RING or NESTED");
        hdu->card("COORDSYS", coordsys(),
                  "Coordinate system, either EQU or GAL");

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Pixel number (0,1,...,num_pixels).
 *
 * HEALPix pixels have all the same solid angle, hence the pix argument is
 * not used.
 ***************************************************************************/
double GHealpix::omega(const int& pix) const
{
    // Return solid angle
    return m_omega;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 *
 * @param[in] pix Sky pixel.
 *
 * HEALPix pixels have all the same solid angle, hence the pix argument is
 * not used.
 ***************************************************************************/
double GHealpix::omega(const GSkyPixel& pix) const
{
    // Return solid angle
    return m_omega;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 ***************************************************************************/
GSkyDir GHealpix::pix2dir(const int& pix) const
{
    // Declare result
    GSkyDir result;
    double  theta = 0.0;
    double  phi   = 0.0;

    // Perform ordering dependent conversion
    switch (m_ordering) {
    case 0:
        pix2ang_ring(pix, &theta, &phi);
        break;
    case 1:
        pix2ang_nest(pix, &theta, &phi);
        break;
    default:
        break;
    }

    // Store coordinate system dependent result
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

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns pixel for a given sky direction
 *
 * @param[in] dir Sky direction.
 ***************************************************************************/
int GHealpix::dir2pix(const GSkyDir& dir) const
{
    // Declare result
    int    pix;
    double z;
    double phi;

    // Compute coordinate system dependent (z,phi)
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
    switch (m_ordering) {
    case 0:
        pix = ang2pix_z_phi_ring(z, phi);
        break;
    case 1:
        pix = ang2pix_z_phi_nest(z, phi);
        break;
    default:
        break;
    }

    // Return pixel index
    return pix;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Sky pixel.
 *
 * This dummy method throws an error when called as sky pixels are not
 * implemented for a HealPix grid. Maybe we should implement them?
 *
 * @todo Think about implementation of sky pixel for HealPix data. We may
 *       use only the first argument, and even carry a usage flag in
 *       GSkyPixel
 ***************************************************************************/
GSkyDir GHealpix::xy2dir(const GSkyPixel& pix) const
{
    // Set error message
    std::string message = "Method not defined for HPX projection.";

    // Throw error
    throw GException::wcs(G_XY2DIR, message);

    // Set dummy return value
    GSkyDir result;
    
    // Return
    return result;
}


/***********************************************************************//**
 * @brief Returns pixel of sky direction
 *
 * @param[in] dir Sky direction.
 *
 * This dummy method throws an error when called as sky pixels are not
 * implemented for a HealPix grid. Maybe we should implement them?
 *
 * @todo Think about implementation of sky pixel for HealPix data. We may
 *       use only the first argument, and even carry a usage flag in
 *       GSkyPixel
 ***************************************************************************/
GSkyPixel GHealpix::dir2xy(const GSkyDir& dir) const
{
    // Set error message
    std::string message = "Method not defined for HPX projection.";

    // Throw error
    throw GException::wcs(G_DIR2XY2, message);

    // Set dummy return value
    GSkyPixel result;
    
    // Return
    return result;
}


/***********************************************************************//**
 * @brief Returns number of pixels
 ***************************************************************************/
int GHealpix::npix(void) const
{
    // Return number of pixels
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns number of divisions of the side of each base pixel.
 ***************************************************************************/
int GHealpix::nside(void) const
{
    // Return nside
    return m_nside;
}


/***********************************************************************//**
 * @brief Returns ordering parameter.
 ***************************************************************************/
std::string GHealpix::ordering(void) const
{
    // Set pixel ordering type
    std::string s_ordering;
    switch (m_ordering) {
    case 0:
        s_ordering = "RING";
        break;
    case 1:
        s_ordering = "NESTED";
        break;
    default:
        s_ordering = "UNKNOWN";
        break;
    }

    // Return ordering
    return s_ordering;
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
        result.append(gammalib::str(m_omega)+" sr");
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
    //m_type        = "HPX";
    m_nside       = 0;
    m_npface      = 0;
    m_ncap        = 0;
    m_order       = 0;
    m_ordering    = 0;
    m_num_pixels  = 0;
    m_fact1       = 0.0;
    m_fact2       = 0.0;
    m_omega       = 0.0;

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
    m_omega      = proj.m_omega;

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
 * @param[in] wcs Pointer to World Coordinate System
 *
 * This method is a helper for the World Coordinate Comparison friends.
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


/***********************************************************************//**
 * @brief Convert nside to order
 *
 * @param[in] nside Number of sides.
 ***************************************************************************/
int GHealpix::nside2order(int nside)
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
    if (ipix < 0 || ipix >= m_num_pixels)
        throw  GException::out_of_range(G_PIX2ANG_RING, ipix, 0, m_num_pixels-1);

    // Handle North Polar cap
    if (ipix < m_ncap) {
        int iring = int(0.5*(1+isqrt(1+2*ipix))); // counted from North pole
        int iphi  = (ipix+1) - 2*iring*(iring-1);
        *theta    = acos(1.0 - (iring*iring) * m_fact2);
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
    if (ipix < 0 || ipix >= m_num_pixels)
        throw GException::out_of_range(G_PIX2ANG_NEST, ipix, 0, m_num_pixels-1);

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
    *theta = acos(z);
    *phi   = (jp - (kshift+1)*0.5) * (gammalib::pihalf / nr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pixels which contains angular coordinates (z,phi)
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
 * @brief Returns pixels which contains angular coordinates (z,phi)
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
 * @brief Integer n that fulfills n*n <= arg < (n+1)*(n+1)
 *
 * @param[in] arg Argument.
 *
 * Returns the integer \a n, which fulfills \a n*n <= arg < (n+1)*(n+1).
 ***************************************************************************/
unsigned int GHealpix::isqrt(unsigned int arg) const
{
    // Return
    return unsigned(std::sqrt(arg+0.5));
}


/***********************************************************************//**
 * @brief Setup of projection
 *
 * Dummy (only implemented for non HPX classes).
 ***************************************************************************/
/*
void GHealpix::prj_set(void)
{
    // Return
    return;
}
*/

/***********************************************************************//**
 * @brief Cartesian-to-spherical deprojection
 *
 * Dummy (only implemented for non HPX classes).
 ***************************************************************************/
/*
int GHealpix::prj_x2s(int nx, int ny, int sxy, int spt, 
                     const double* x, const double* y,
                     double* phi, double* theta, int* stat) const
{
    // Return
    return 0;
}
*/

/***********************************************************************//**
 * @brief Generic spherical-to-Cartesian projection
 *
 * Dummy (only implemented for non HPX classes).
 ***************************************************************************/
/*
int GHealpix::prj_s2x(int nphi, int ntheta, int spt, int sxy,
                     const double* phi, const double* theta,
                     double* x, double* y, int* stat) const
{
    // Return
    return 0;
}
*/
