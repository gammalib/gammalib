/***************************************************************************
 *             GHealpix.cpp  -  Healpix sky representation class           *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
 * @file GHealpix.cpp
 * @brief GHealpix class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <math.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GHealpix.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_LOAD         "GHealpix::load(const GFitsHDU*)"
#define G_PIX2ANG_RING "GHealpix::pix2ang_ring(int,double*,double*)"
#define G_PIX2ANG_NEST "GHealpix::pix2ang_nest(int,double*,double*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const int jrll[12] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
const int jpll[12] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};

/* __ Static conversion arrays ___________________________________________ */
static int pix2x[1024];
static int pix2y[1024];


/*==========================================================================
 =                                                                         =
 =                     GHealpix constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GHealpix::GHealpix()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param pixels GHealpix instance which should be used for construction
 ***************************************************************************/
GHealpix::GHealpix(const GHealpix& pixels)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(pixels);

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
 =                          GHealpix operators                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] pixels GHealpix instance to be assigned
 ***************************************************************************/
GHealpix& GHealpix::operator= (const GHealpix& pixels)
{
    // Execute only if object is not identical
    if (this != &pixels) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(pixels);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GHealpix public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Load Healpix data from table.
 *
 * @param[in] hdu FITS HDU containing the Healpix data
 *
 * @exception GException::healpix
 *            Unable to load Healpix data from table.
 ***************************************************************************/
void GHealpix::load(const GFitsHDU* hdu)
{
    // Free memory and initialise members
    free_members();
    init_members();
    
    // Check if we have a healpix representation
    if (hdu->card("PIXTYPE")->string() != "HEALPIX")
        throw GException::healpix(G_LOAD, "HDU does not contain Healpix data");
    
    // Get Healpix resolution and determine number of pixels and solid angle
    m_nside      = hdu->card("NSIDE")->integer();
    m_num_pixels = 12 * m_nside * m_nside;
    m_omega      = fourpi / m_num_pixels;

    // Get ordering scheme from ORDERING keyword
    std::string ordering;
    try {
        ordering = hdu->card("ORDERING")->string();
    }
    catch (GException::fits_key_not_found &e) {
        m_order = -1; // Flag unknown ordering
    }
    
    // Decode ordering scheme string
    if (ordering == "RING")
        m_order = 0;
    else if (ordering == "NESTED")
        m_order = 1;
    else
        throw GException::healpix(G_LOAD, "Invalid ordering");

    // Get coordinate system. First search for the COORDSYS keyword.
    // If not found then search for the HIER_CRD keyword. This has
    // been used in older versions of LAT exposure cubes ...
    std::string coordsys;
    try {
        coordsys = hdu->card("COORDSYS")->string();
    }
    catch (GException::fits_key_not_found &e) {
        try {
            coordsys = hdu->card("HIER_CRD")->string();
        }
        catch (GException::fits_key_not_found &e) {
            m_coordsys = -1; // Flag coordinate system as unknown
        }
    }
    
    // Decode coordinate system string
    if (coordsys.find("EQU") == 0)
        m_coordsys = 0;
    else if (coordsys.find("GAL") == 0)
        m_coordsys = 1;
    else
        throw GException::healpix(G_LOAD, "Invalid coordinate system");
    
    // Continue only of we have pixels
    if (m_num_pixels > 0) {
    
        // Get first column
        GFitsTableCol* col = hdu->column(0);

        // Check column consistency
        if (col->length() != m_num_pixels)
            throw GException::healpix(G_LOAD, 
                  "NSIDE inconsistent with number of rows");

        // Extract vector size of each pixel
        m_size_pixels = col->number();

        // If there are pixels then load them
        int size = m_num_pixels * m_size_pixels;
        if (size > 0) {
            m_pixels = new double[size];
            double* ptr = m_pixels;
            for (int row = 0; row < m_num_pixels; ++row) {
                for (int inx = 0; inx < m_size_pixels; ++inx, ++ptr)
                    *ptr = col->real(row,inx);
            }
        } // endif: there were pixels to load

        // Compute sky direction for each pixel
        m_dir = new GSkyDir[m_num_pixels];
        for (int i = 0; i < m_num_pixels; ++i)
            m_dir[i] = pix2ang(i);
    
    } // endif: we had pixels
        
    // Return
    return;
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
 * @brief Returns number of pixels.
 ***************************************************************************/
int GHealpix::num_pixels(void) const
{
    // Return nside
    return m_num_pixels;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel
 ***************************************************************************/
double GHealpix::omega(void) const
{
    // Return solid angle
    return m_omega;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] ipix Pixel number (0,1,...,num_pixels)
 ***************************************************************************/
GSkyDir GHealpix::pix2ang(const int& ipix)
{
    // Declare result
    GSkyDir result;
    double  theta = 0.0;
    double  phi   = 0.0;
    
    // Perform ordering dependent conversion
    switch (m_order) {
    case 0:
        pix2ang_ring(ipix, &theta, &phi);
        break;
    case 1:
        pix2ang_nest(ipix, &theta, &phi);
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


/*==========================================================================
 =                                                                         =
 =                         GHealpix private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GHealpix::init_members(void)
{
    // Initialise members
    m_nside       = 0;
    m_order       = 0;
    m_coordsys    = 0;
    m_num_pixels  = 0;
    m_size_pixels = 0;
    m_omega       = 0.0;
    m_pixels      = NULL;
    m_dir         = NULL;
    
    // Initiates the array for the pixel number -> (x,y) mapping
    mk_pix2xy();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pixels GHealpix instance from which members should be copied
 ***************************************************************************/
void GHealpix::copy_members(const GHealpix& pixels)
{
    // Copy attributes
    m_nside       = pixels.m_nside;
    m_order       = pixels.m_order;
    m_coordsys    = pixels.m_coordsys;
    m_num_pixels  = pixels.m_num_pixels;
    m_size_pixels = pixels.m_size_pixels;
    m_omega       = pixels.m_omega;
    
    // Copy arrays
    if (m_num_pixels > 0) {
        if (m_size_pixels > 0 && pixels.m_pixels != NULL) {
            int size = m_num_pixels*m_size_pixels;
            m_pixels = new double[size];
            memcpy(m_pixels, pixels.m_pixels, size*sizeof(double));
        }
        if (pixels.m_dir != NULL) {
            m_dir = new GSkyDir[m_num_pixels];
            for (int i = 0; i < m_num_pixels; ++i)
                m_dir[i] = pixels.m_dir[i];
        }
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GHealpix::free_members(void)
{
    // Free memory
    if (m_pixels != NULL) delete [] m_pixels;
    if (m_dir    != NULL) delete [] m_dir;

    // Mark memory as free
    m_pixels = NULL;
    m_dir    = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Healpix representation
 ***************************************************************************/
GHealpix* GHealpix::clone(void) const
{
    return new GHealpix(*this);
}


/***********************************************************************//**
 * @brief Constructs pixel to (x,y) conversion array
 *
 * Constructs the array giving x and y in the face from pixel number for the
 * nested (quad-cube like) ordering of pixels. The bits corresponding to x 
 * and y are interleaved in the pixel number one breaks up the pixel number
 * by even and odd bits.
 ***************************************************************************/
void GHealpix::mk_pix2xy(void)
{
    // Loop over all pixels of array
    for (int kpix = 0; kpix < 1024; ++kpix) {
    
        // Setup
        int jpix = kpix;
        int IX   = 0;
        int IY   = 0;
        int IP   = 1;                    // Bit position (in x and y)
        
        // Go through all the bits
        while (jpix != 0) {
            int ID = (int)fmod(jpix,2);  //  Bit value (in kpix), goes in ix
            jpix   = jpix/2;
            IX     = ID*IP+IX;
            ID     = (int)fmod(jpix,2);  //  Bit value (in kpix), goes in iy
            jpix   = jpix/2;
            IY     = ID*IP+IY;
            IP     = 2*IP;               //  Next bit (in x and y)
        }

        // Set array elements
        pix2x[kpix] = IX;                //  in {0,31}
        pix2y[kpix] = IY;                //  in {0,31}
        
    }
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert pixel index to (theta,phi) angles for ring ordering
 *
 * @param[in] ipix Pixel index for which (theta,phi) are to be computed.
 * @param[out] theta Pointer to result zenith angle in radians.
 * @param[out] phi Pointer to result azimuth angle in radians.
 *
 * @exception GException::out_of_range Pixel index is out of range.
 ***************************************************************************/
void GHealpix::pix2ang_ring(int ipix, double* theta, double* phi)
{
    // Check if ipix is in range
    if (ipix < 0 || ipix >= m_num_pixels)
        throw  GException::out_of_range(G_PIX2ANG_RING, ipix, 0, m_num_pixels-1);

    // Setup
    int ipix1 = ipix + 1;     // in {1, m_num_pixels}
    int nl2   = 2 * m_nside;
    int nl4   = 4 * m_nside;
    
    // Determine points in each polar cap (0 for nside=1)
    int ncap  = 2 * m_nside * (m_nside-1);
    
    double fact1 = 1.5 * m_nside;
    double fact2 = 3.0 * m_nside * m_nside;
  
    // Handle North Polar cap
    if (ipix1 <= ncap) {
        double hip   = ipix1/2.;
        double fihip = floor(hip);
        int    iring = (int)floor(sqrt(hip - sqrt(fihip)))+1; 
        int    iphi  = ipix1 - 2*iring*(iring - 1);
        *theta       = acos(1. - iring*iring / fact2);
        *phi         = (1.*iphi - 0.5) * pi/(2.*iring);
    }
    
    // Handle Equatorial region
    else if(ipix1 <= nl2*(5*m_nside+1)) {
        int ip       = ipix1 - ncap - 1;
        int iring    = (int)floor(ip/nl4) + m_nside;
        int iphi     = (int)fmod(ip,nl4) + 1;
        double fodd  = 0.5 * (1 + fmod((double)(iring+m_nside),2));
        *theta       = acos((nl2 - iring) / fact1);
        *phi         = (1.*iphi - fodd) * pi/(2.*m_nside);
    }
    
    // Handle South Polar cap
    else {
        int    ip    = m_num_pixels - ipix1 + 1;
        double hip   = ip/2.;
        double fihip = floor(hip);
        int    iring = (int)floor(sqrt( hip - sqrt(fihip))) + 1;
        int    iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));
        *theta       = acos( -1. + iring*iring / fact2 );
        *phi         = (1.*iphi - 0.5) * pi/(2.*iring);
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
 * @exception GException::out_of_range Pixel index is out of range.
 ***************************************************************************/
void GHealpix::pix2ang_nest(int ipix, double* theta, double* phi)
{
    // Check if ipix is in range
    if (ipix < 0 || ipix >= m_num_pixels)
        throw GException::out_of_range(G_PIX2ANG_NEST, ipix, 0, m_num_pixels-1);

    // ...
    double fn    = 1. * m_nside;
    double fact1 = 1. / (3.*fn*fn);
    double fact2 = 2. / (3.*fn);
    int    nl4   = 4 * m_nside;

    // Finds the face, and the number in the face
    int npface   = m_nside*m_nside;
    int face_num = ipix/npface;            // Face number in {0,11}
    int ipf      = (int)fmod(ipix,npface); // Pixel number in the face {0,npface-1}

    // Finds the x,y on the face from the pixel number
    // (starting from the lowest corner)
    int ip_low   = (int)fmod(ipf,1024);      // Content of the last 10 bits
    int ip_trunc = ipf/1024;                 // Truncation of the last 10 bits
    int ip_med   = (int)fmod(ip_trunc,1024); // Content of the next 10 bits
    int ip_hi    = ip_trunc/1024;            // Content of the high weight 10 bits
    int  ix      = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
    int  iy      = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

    // Transforms this in (horizontal, vertical) coordinates
    int jrt = ix + iy;   // 'vertical' in {0,2*(m_nside-1)}
    int jpt = ix - iy;   // 'horizontal' in {-m_nside+1,m_nside-1}

    // Computes the z coordinate on the sphere
    int    jr     = jrll[face_num]*m_nside - jrt - 1;
    int    nr     = m_nside;             // Equatorial region (the most frequent)
    double z      = (2*m_nside-jr)*fact2;
    int    kshift = (int)fmod(jr - m_nside, 2);
    
    // North pole region
    if (jr < m_nside) {
        nr     = jr;
        z      = 1. - nr*nr*fact1;
        kshift = 0;
    }
    
    // South pole region
    else if (jr > 3*m_nside) {
        nr     = nl4 - jr;
        z      = - 1. + nr*nr*fact1;
        kshift = 0;
    }
    
    // Computes theta
    *theta = acos(z);
      
    // Computes the phi coordinate on the sphere, in [0,2Pi]
    int jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
    if (jp > nl4) jp = jp - nl4;
    if (jp <   1) jp = jp + nl4;
    
    // Computes Phi
    *phi = (jp - (kshift+1)*0.5) * (pihalf / nr);
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GHealpix friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] column Healpix array to put in output stream
 ***************************************************************************/
ostream& operator<< (ostream& os, const GHealpix& pixels)
{
    // Put header in stream
    os << "=== GHealpix ===" << endl;
    os << " Number of base pixel div. .: " << pixels.m_nside << endl;
    os << " Number of pixels ..........: " << pixels.m_num_pixels << endl;
    os << " Pixel vector size .........: " << pixels.m_size_pixels << endl;
    os << " Solid angle ...............: " << scientific << pixels.m_omega 
       << fixed << endl;
    
    // Put ordering in stream
    os << " Ordering ..................: ";
    switch (pixels.m_order) {
    case 0:
        os << "Ring" << endl;
        break;
    case 1:
        os << "Nested" << endl;
        break;
    case -1:
        os << "*** Unknown ***" << endl;
        break;
    default:
        os << "*** Invalid ***" << endl;
        break;
    }

    // Put coordinate system in stream
    os << " Coordinate system .........: ";
    switch (pixels.m_coordsys) {
    case 0:
        os << "Equatorial (RA,Dec)" << endl;
        break;
    case 1:
        os << "Galactic (l,b)" << endl;
        break;
    case -1:
        os << "*** Unknown ***" << endl;
        break;
    default:
        os << "*** Invalid ***" << endl;
        break;
    }

    for (int i = 0; i < pixels.m_num_pixels; ++i)
        os << pixels.m_dir[i] << " ";

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                     Other functions used by GHealpix                    =
 =                                                                         =
 ==========================================================================*/
