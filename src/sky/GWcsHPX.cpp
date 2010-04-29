/***************************************************************************
 *                 GWcsHPX.cpp  -  Healpix projection class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GWcsHPX.cpp
 * @brief GWcsHPX class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcsHPX.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCT           "GWcsHPX::GWcsHPX(int,std::string,std::string)"
#define G_READ                                     "GWcsHPX::read(GFitsHDU*)"
#define G_XY2DIR                                 "GWcsHPX::xy2dir(GSkyPixel)"
#define G_DIR2XY                                   "GWcsHPX::dir2xy(GSkyDir)"
#define G_PIX2ANG_RING           "GWcsHPX::pix2ang_ring(int,double*,double*)"
#define G_PIX2ANG_NEST           "GWcsHPX::pix2ang_nest(int,double*,double*)"
#define G_ORDERING_SET                       "GWcsHPX::coordsys(std::string)"

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


/*==========================================================================
 =                                                                         =
 =                      GWcsHPX constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcsHPX::GWcsHPX(void) : GWcs()
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
 * @param[in] dimension Vector dimension of pixels.
 *
 * @exception GException::wcs_hpx_bad_nside 
 *            Invalid nside parameter.
 * @exception GException::wcs_bad_coords 
 *            Invalid coordsys parameter.
 * @exception GException::wcs_hpx_bad_ordering 
 *            Invalid ordering parameter.
 ***************************************************************************/
GWcsHPX::GWcsHPX(const int& nside, const std::string& order,
                 const std::string& coords) : GWcs()
{
    // Initialise class members
    init_members();

    // Check nside parameter (power of 2 between 1 and 8192)
    if (nside2order(nside) == -1)
        throw GException::wcs_hpx_bad_nside(G_CONSTRUCT, nside);

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
    m_omega      = fourpi / m_num_pixels;
    m_order      = nside2order(m_nside);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor from FITS HDU table
 *
 * @param[in] hdu Pointer to FITS HDU.
 ***************************************************************************/
GWcsHPX::GWcsHPX(const GFitsHDU* hdu) : GWcs()
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
 * @param[in] wcs GWcsHPX instance which should be used for construction.
 ***************************************************************************/
GWcsHPX::GWcsHPX(const GWcsHPX& wcs) : GWcs(wcs)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(wcs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GWcsHPX::~GWcsHPX()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            GWcsHPX operators                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] wcs GWcsHPX instance to be assigned.
 ***************************************************************************/
GWcsHPX& GWcsHPX::operator= (const GWcsHPX& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

        // Copy base class members
        this->GWcs::operator=(wcs);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(wcs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GWcsHPX public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns WCS type
 ***************************************************************************/
std::string GWcsHPX::type(void) const
{
    // Return Healix type
    return "HPX";
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
void GWcsHPX::read(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();

    // Check if we have a healpix representation
    if (hdu->card("PIXTYPE")->string() != "HEALPIX")
        throw GException::wcs(G_READ, "HDU does not contain Healpix data");

    // Get pixel ordering
    std::string order = hdu->card("ORDERING")->string();

    // Get coordinate system.
    // First search for HIER_CRD keyword (this has been used in older
    // versions of LAT exposure cubes). If not found then search for standard
    // COORDSYS keyword.
    std::string coords;
    try {
        coords = hdu->card("HIER_CRD")->string();
    }
    catch (GException::fits_key_not_found &e) {
        coords = hdu->card("COORDSYS")->string();
    }

    // Set coordinate system
    coordsys(coords);

    // Set pixel ordering
    ordering(order);

    // Get Healpix resolution and determine number of pixels and solid angle
    m_nside      = hdu->card("NSIDE")->integer();
    m_npface     = m_nside * m_nside;
    m_ncap       = 2 * (m_npface - m_nside);
    m_num_pixels = 12 * m_npface;
    m_fact2      = 4.0 / m_num_pixels;
    m_fact1      = 2 * m_nside * m_fact2;
    m_omega      = fourpi / m_num_pixels;
    m_order      = nside2order(m_nside);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write Healpix definiton into FITS HDU
 *
 * @param[in] hdu FITS HDU to which the Healpix definition will be written.
 *
 * Writes or updates the following keywords in the FITS HDU:
 * EXTNAME  = HEALPIX
 * PIXTYPE  = HEALPIX
 * NSIDE    = nside() = m_nside
 * FIRSTPIX = 0
 * LASTPIX  = npix()-1 = m_num_pixels-1
 * ORDERING = ordering()
 * COORDSYS = coordsys()
 ***************************************************************************/
void GWcsHPX::write(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
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
double GWcsHPX::omega(const int& pix) const
{
    // Return solid angle
    return m_omega;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 ***************************************************************************/
GSkyDir GWcsHPX::pix2dir(const int& pix)
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
        result.radec(phi, pihalf-theta);
        break;
    case 1:
        result.lb(phi, pihalf-theta);
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
int GWcsHPX::dir2pix(GSkyDir dir) const
{
    // Declare result
    int    pix;
    double z;
    double phi;

    // Compute coordinate system dependent (z,phi)
    switch (m_coordsys) {
    case 0:
        z   = cos(pihalf-dir.dec());
        phi = dir.ra();
        break;
    case 1:
        z   = cos(pihalf-dir.b());
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
 * @brief Returns number of pixels
 ***************************************************************************/
int GWcsHPX::npix(void) const
{
    // Return number of pixels
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns number of divisions of the side of each base pixel.
 ***************************************************************************/
int GWcsHPX::nside(void) const
{
    // Return nside
    return m_nside;
}


/***********************************************************************//**
 * @brief Returns ordering parameter.
 ***************************************************************************/
std::string GWcsHPX::ordering(void) const
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
void GWcsHPX::ordering(const std::string& ordering)
{
    // Convert argument to upper case
    std::string uordering = toupper(ordering);

    // Set pixel ordering
    if (uordering == "RING")
        m_ordering = 0;
    else if (uordering == "NESTED" || uordering == "NEST")
        m_ordering = 1;
    else
        throw GException::wcs_hpx_bad_ordering(G_ORDERING_SET, ordering);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GWcsHPX private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GWcsHPX::init_members(void)
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
 * @param[in] wcs GWcsHPX instance from which members should be copied.
 ***************************************************************************/
void GWcsHPX::copy_members(const GWcsHPX& wcs)
{
    // Copy attributes
    m_nside       = wcs.m_nside;
    m_npface      = wcs.m_npface;
    m_ncap        = wcs.m_ncap;
    m_order       = wcs.m_order;
    m_ordering    = wcs.m_ordering;
    m_num_pixels  = wcs.m_num_pixels;
    m_fact1       = wcs.m_fact1;
    m_fact2       = wcs.m_fact2;
    m_omega       = wcs.m_omega;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsHPX::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GWcsHPX* GWcsHPX::clone(void) const
{
    return new GWcsHPX(*this);
}


/***********************************************************************//**
 * @brief Convert nside to order
 *
 * @param[in] nside Number of sides.
 ***************************************************************************/
int GWcsHPX::nside2order(int nside)
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
void GWcsHPX::pix2xy(const int& ipix, int* x, int* y)
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
int GWcsHPX::xy2pix(int x, int y) const
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
void GWcsHPX::pix2ang_ring(int ipix, double* theta, double* phi)
{
    // Check if ipix is in range
    if (ipix < 0 || ipix >= m_num_pixels)
        throw  GException::out_of_range(G_PIX2ANG_RING, ipix, 0, m_num_pixels-1);

    // Handle North Polar cap
    if (ipix < m_ncap) {
        int iring = int(0.5*(1+isqrt(1+2*ipix))); // counted from North pole
        int iphi  = (ipix+1) - 2*iring*(iring-1);
        *theta    = acos(1.0 - (iring*iring) * m_fact2);
        *phi      = (iphi - 0.5) * pi/(2.0*iring);
    }

    // Handle Equatorial region
    else if (ipix < (m_num_pixels - m_ncap)) {
        int    ip    = ipix - m_ncap;
        int    iring = ip/(4*m_nside) + m_nside;   // counted from North pole
        int    iphi  = ip%(4*m_nside) + 1;
        double fodd  = ((iring+m_nside)&1) ? 1 : 0.5;
        int    nl2   = 2*m_nside;
        *theta       = acos((nl2 - iring) * m_fact1);
        *phi         = (iphi - fodd) * pi/nl2;
    }

    // Handle South Polar cap
    else {
        int ip    = m_num_pixels - ipix;
        int iring = int(0.5*(1+isqrt(2*ip-1)));    // Counted from South pole
        int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
        *theta    = acos(-1.0 + (iring*iring) * m_fact2);
        *phi      = (iphi - 0.5) * pi/(2.*iring);
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
void GWcsHPX::pix2ang_nest(int ipix, double* theta, double* phi)
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
    *phi   = (jp - (kshift+1)*0.5) * (pihalf / nr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pixels which contains angular coordinates (z,phi)
 *
 * @param[in] z Cosine of zenith angle - cos(theta).
 * @param[in] phi Azimuth angle in radians.
 ***************************************************************************/
int GWcsHPX::ang2pix_z_phi_ring(double z, double phi) const
{
    // Initialise pixel
    int ipix = 0;

    // Setup
    double za = fabs(z);
    double tt = modulo(phi,twopi) * inv_pihalf; // in [0,4)

    // Equatorial region
    if (za <= twothird) {
        double temp1  = m_nside*(0.5+tt);
        double temp2  = m_nside*z*0.75;
        int    jp     = int(temp1-temp2);           // index of ascending edge line
        int    jm     = int(temp1+temp2);           // index of descending edge line
        int    ir     = m_nside + 1 + jp - jm;      // in {1,2n+1}
        int    kshift = 1 - (ir & 1);               // kshift=1 if ir even, 0 otherwise
        int    ip     = (jp+jm-m_nside+kshift+1)/2; // in {0,4n-1}
        ip            = int(modulo(ip,4*m_nside));
        ipix          = m_ncap + (ir-1)*4*m_nside + ip;
    }

    // North & South polar caps
    else {
        double tp  = tt - int(tt);
        double tmp = m_nside * sqrt(3*(1-za));
        int    jp  = int(tp*tmp);       // increasing edge line index
        int    jm  = int((1.0-tp)*tmp); // decreasing edge line index
        int    ir  = jp + jm + 1;       // ring number counted from the closest pole
        int    ip  = int(tt*ir);        // in {0,4*ir-1}
        ip = int(modulo(ip,4*ir));
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
int GWcsHPX::ang2pix_z_phi_nest(double z, double phi) const
{
    // Initialise face and pixel numbers
    int face_num;
    int ix;
    int iy;

    // Setup
    double za = fabs(z);
    double tt = modulo(phi,twopi) * inv_pihalf; // in [0,4)

    // Equatorial region
    if (za <= twothird) {
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
        double tmp = ns_max*sqrt(3*(1-za));
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
unsigned int GWcsHPX::isqrt(unsigned int arg)
{
    // Return
    return unsigned(sqrt(arg+0.5));
}


/*==========================================================================
 =                                                                         =
 =                              GWcsHPX friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] wcs Healpix WCS definition to put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GWcsHPX& wcs)
{
    // Put header in stream
    os << "=== GWcsHPX ===" << std::endl;
    os << " Nside (number of divisions): " << wcs.m_nside << std::endl;
    os << " Npface (pixels per face) ..: " << wcs.m_npface << std::endl;
    os << " Ncap (number of cap pixels): " << wcs.m_ncap << std::endl;
    os << " Npix (number of pixels) ...: " << wcs.m_num_pixels << std::endl;
    os << " Order .....................: " << wcs.m_order << std::endl;
    os << " Solid angle ...............: " << std::scientific << wcs.m_omega
       << std::fixed << " sr" << std::endl;
    os << " Ordering ..................: " << wcs.ordering() << std::endl;
    os << " Coordinate system .........: " << wcs.coordsys();

    // Return output stream
    return os;
}
