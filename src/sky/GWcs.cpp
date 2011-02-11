/***************************************************************************
 *           GWcs.cpp  -  World Coordinate System virtual base class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GWcs.cpp
 * @brief World Coordinate System virtual base class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COORDSYS_SET                          "GWcs::coordsys(std::string)"
#define G_OMEGA1                                           "GWcs::omega(int)"
#define G_OMEGA2                                     "GWcs::omega(GSkyPixel)"
#define G_XY2DIR                                    "GWcs::xy2dir(GSkyPixel)"
#define G_DIR2XY                                      "GWcs::dir2xy(GSkyDir)"
#define G_WCS_SET      "GWcs::GWcs(std::string,double,double,double,double,"\
                                                             "double,double)"
#define G_WCS_READ                                "GWcs::wcs_read(GFitsHDU*)"
#define G_WCS_WRITE                              "GWcs::wcs_write(GFitsHDU*)"
#define G_WCS_CRVAL1                                     "GWcs::wcs_crval1()"
#define G_WCS_CRVAL2                                     "GWcs::wcs_crval1()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DIR2XY_DEBUG                                      // Debug dir2xy
//#define G_XY2DIR_DEBUG                                      // Debug xy2dir

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Static conversion arrays ___________________________________________ */

/*==========================================================================
 =                                                                         =
 =                       GWcs constructors/destructors                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcs::GWcs(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Standard WCS sky map constructor
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel (first pixel is 1).
 * @param[in] crpix2 Y index of reference pixel (first pixel is 1).
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 *
 * Construct standard WCS sky map from standard definition parameters. This
 * method
 ***************************************************************************/
GWcs::GWcs(const std::string& coords,
           const double& crval1, const double& crval2,
           const double& crpix1, const double& crpix2,
           const double& cdelt1, const double& cdelt2)

{
    // Initialise class members
    init_members();

    // Initialise reverse flag
    m_reverse = 0;

    // Set standard parameters
    wcs_set(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param wcs World Coordinate System.
 ***************************************************************************/
GWcs::GWcs(const GWcs& wcs)
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
GWcs::~GWcs(void)
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
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
GWcs& GWcs::operator= (const GWcs& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return WCS type
 ***************************************************************************/
std::string GWcs::type(void) const
{
    // Return Healix type
    return m_type;
}


/***********************************************************************//**
 * @brief Returns coordinate system.
 *
 * Returns either 'EQU' or 'GAL'.
 ***************************************************************************/
std::string GWcs::coordsys(void) const
{
    // Set coordinate system
    std::string s_coordsys;
    switch (m_coordsys) {
    case 0:
        s_coordsys = "EQU";
        break;
    case 1:
        s_coordsys = "GAL";
        break;
    default:
        s_coordsys = "UNKNOWN";
        break;
    }

    // Return coordinate system
    return s_coordsys;
}


/***********************************************************************//**
 * @brief Set coordinate system
 *
 * @param[in] coordsys Coordinate system (EQU/CEL/E/C or GAL/G)
 *
 * @exception GException::wcs_bad_coords
 *            Invalid coordsys parameter.
 *
 * Set coordinate system from std::string. Each of the following is 
 * interpreted as celestial coordinate system: EQU, CEL, E, C. Each of the
 * following is interpreted as galactic coordinate system: GAL, G.
 ***************************************************************************/
void GWcs::coordsys(const std::string& coordsys)
{
    // Convert argument to upper case
    std::string ucoordsys = toupper(coordsys);

    // Set coordinate system
    if (ucoordsys == "EQU" || ucoordsys == "CEL" || ucoordsys == "E" ||
        ucoordsys == "C")
        m_coordsys = 0;
    else if (ucoordsys == "GAL" || ucoordsys == "G")
        m_coordsys = 1;
    else
        throw GException::wcs_bad_coords(G_COORDSYS_SET, coordsys);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel (DUMMY METHOD)
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 *
 * @exception GException::wcs
 *            Method not defined for projection.
 *
 * This method should be overloaded by the appropriate method of the derived
 * class. It thus should never get called.
 ***************************************************************************/
double GWcs::omega(const int& pix) const
{
    // Set error message
    std::string message = "GWcs method not defined for " +
                          this->type() + " projection.";

    // Throw error
    throw GException::wcs(G_OMEGA1, message);

    // Define constant direction
    const double omega = 0.0;

    // Return solid angle
    return omega;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel (DUMMY METHOD)
 *
 * @param[in] pix Sky pixel.
 *
 * @exception GException::wcs
 *            Method not defined for projection.
 *
 * This method should be overloaded by the appropriate method of the derived
 * class. It thus should never get called.
 ***************************************************************************/
double GWcs::omega(const GSkyPixel& pix) const
{
    // Set error message
    std::string message = "GWcs method not defined for " +
                          this->type() + " projection.";

    // Throw error
    throw GException::wcs(G_OMEGA2, message);

    // Define constant direction
    const double omega = 0.0;

    // Return solid angle
    return omega;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel (DUMMY METHOD)
 *
 * @param[in] pix Pixel number (0,1,...,m_num_pixels).
 *
 * @exception GException::wcs
 *            Method not defined for projection.
 *
 * This method should be overloaded by the appropriate method of the derived
 * class. It thus should never get called.
 ***************************************************************************/
GSkyDir GWcs::pix2dir(const int& pix) const
{
    // Set error message
    std::string message = "GWcs method not defined for " +
                          this->type() + " projection.";

    // Throw error
    throw GException::wcs(G_XY2DIR, message);

    // Define constant direction
    const GSkyDir dir;

    // Return
    return dir;
}


/***********************************************************************//**
 * @brief Returns pixel for sky direction (DUMMY METHOD)
 *
 * @param[in] dir Sky direction.
 *
 * @exception GException::wcs
 *            Method not defined for projection.
 *
 * This method should be overloaded by the appropriate method of the derived
 * class. It thus should never get called.
 ***************************************************************************/
int GWcs::dir2pix(GSkyDir dir) const
{
    // Set error message
    std::string message = "GWcs method not defined for " +
                          this->type() + " projection.";

    // Throw error
    throw GException::wcs(G_DIR2XY, message);

    // Return
    return 0;
}


/***********************************************************************//**
 * @brief Returns sky direction of pixel
 *
 * @param[in] pix Sky pixel.
 *
 * Note that pixel indices start from 0.
 ***************************************************************************/
GSkyDir GWcs::xy2dir(const GSkyPixel& pix) const
{
    // Set constants
    const int __permutation[2] = {1,0};

    // Determine offset
    GVector offset = GVector(pix.x(), pix.y());

    // Determine xy
    GVector xy = m_cd * (offset - m_refpix);

    // Swap result if coordinates are reversed
    if (m_reverse)
        xy = perm(xy, __permutation);

    // Perform CAR map projection
    GVector native = xy;
    nat2std(&native);

    // Get sky direction for native coordinates
    GSkyDir dir = wcs_native2dir(native);

    // Debug: Dump transformation steps
    #if defined(G_XY2DIR_DEBUG)
    std::cout << "xy2dir: pixel=" << offset << " xy=" << xy
              << " native=" << native << " dir=" << dir << std::endl;
    #endif

    // Return
    return dir;
}


/***********************************************************************//**
 * @brief Returns pixel of sky direction
 *
 * @param[in] dir Sky direction.
 *
 * Note that pixel indices start from 0.
 ***************************************************************************/
GSkyPixel GWcs::dir2xy(GSkyDir dir) const
{
    // Set constants
    const int __permutation[2] = {1,0};

    // Get native coordinates for sky direction in radians
    GVector native = wcs_dir2native(dir);

    // Perform CAR map projection
    GVector xy = native;
    std2nat(&xy);

    // Swap result if coordinates are reversed
    if (m_reverse)
        xy = perm(xy, __permutation);

    // Determine pixel offset
    GVector offset = m_invcd * xy + m_refpix;

    // Set sky pixel
    GSkyPixel pixel(offset[0], offset[1]);

    // Debug: Dump transformation steps
    #if defined(G_DIR2XY_DEBUG)
    std::cout << "dir2xy: dir=" << dir << " native=" << native
              << " xy=" << xy << " pixel=" << offset << std::endl;
    #endif

    // Return sky pixel
    return pixel;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Returns native coordinates (in radians) for sky direction
 *
 * @param[in] dir Sky direction.
 ***************************************************************************/
GVector GWcs::wcs_dir2native(GSkyDir dir) const
{
    // Get sky direction in radians and map coordinate system
    double phi1;
    double theta1;
    if (m_coordsys == 0) {
        phi1   = dir.ra();
        theta1 = dir.dec();
    }
    else {
        phi1   = dir.l();
        theta1 = dir.b();
    }

    // Define right hand side vector for matrix equation
    double  ct = cos(theta1);
    double  x0 = ct * cos(phi1);
    double  x1 = ct * sin(phi1);
    double  x2 = sin(theta1);
    GVector x(x0,x1,x2);

    // Find solution
    GVector b = m_rot * x;

    // Account for possible roundoff in 3rd element
    if (b[2] < -1.0) b[2] = -1.0;
    if (b[2] > +1.0) b[2] = +1.0;

    // Compute native coordinates in radians
    double  theta = asin(b[2]);
    double  phi   = atan2(b[1],b[0]);
    GVector native(phi, theta);

    // Return native coordinates
    return native;
}


/***********************************************************************//**
 * @brief Returns sky direction for native coordinates
 *
 * @param[in] native Native coordinates (in radians).
 ***************************************************************************/
GSkyDir GWcs::wcs_native2dir(GVector native) const
{
    // Get native coordinates
    double phi1   = native[0];
    double theta1 = native[1];

    // Define right hand side vector for matrix equation
    double  ct = cos(theta1);
    double  x0 = ct * cos(phi1);
    double  x1 = ct * sin(phi1);
    double  x2 = sin(theta1);
    GVector x(x0,x1,x2);

    // Find solution
    GVector b = m_trot * x;

    // Account for possible roundoff in 3rd element
    if (b[2] < -1.0) b[2] = -1.0;
    if (b[2] > +1.0) b[2] = +1.0;

    // Compute sky coordinates in radians
    double  theta = asin(b[2]);
    double  phi   = atan2(b[1],b[0]);

    // Assign sky direction
    GSkyDir dir;
    if (m_coordsys == 0)
        dir.radec(phi, theta);
    else
        dir.lb(phi, theta);

    // Return sky direction
    return dir;
}


/***********************************************************************//**
 * @brief Initialises derived WCS parameters
 *
 * @param[in] theta0 Native latitude of the fiducial point.
 ***************************************************************************/
void GWcs::wcs_init(const double& theta0)
{
    // Set native latitude of the fiducial point
    m_theta0 = theta0;

    // Compute the coordinates of native pole for non-polar projection
    m_native_pole = wcs_getpole(m_theta0);

    // Get matrix for rotation between sky and native coordinates
    m_rot = wcs_get_rot();

    // Get transpose matrix for inverse rotation
    m_trot = transpose(m_rot);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute the coordinates of native pole for non-polar projection
 *
 * @param[in] theta0 Native latitude of the fiducial point.
 *
 * For non-polar (cylindrical or conic) projections, the native pole is
 * not at the reference point, and this method is used to determine the
 * position of the native pole. See section 2.4 of the paper
 * "Representation of Celestial Coordinates in FITS" by 
 * Calabretta & Greisen (2002) A&A, 395, 1077, also available at
 * http://www.aoc.nrao.edu/~egreisen.
 * This method returns the coordinates of the native pole for non-polar
 * projections in spherical coordinates. Units are radians.
 * The code has largely be inspired by the IDL routine wcs_getpole.pro
 * (update of May 1998).
 ***************************************************************************/
GVector GWcs::wcs_getpole(const double& theta0)
{
    // Declare result
    double alpha_p;
    double delta_p;

    // Compute crval in radians
    double alpha_0 = m_refval[0] * deg2rad;
    double delta_0 = m_refval[1] * deg2rad;

    // Get longitude and latitude of pole
    double lonpole = m_npole[0];
    double latpole = m_npole[1];

    // If theta0=90 then the coordinate of the native pole are given by crval
    if (theta0 == 90.0) {
        alpha_p = alpha_0;
        delta_p = delta_0;
    }

    // ... otherwise perform computations
    else {
        // Compute some angles
        double phi_p = lonpole*deg2rad;
        double sp    = sin(phi_p);
        double cp    = cos(phi_p);
        double sd    = sin(delta_0);
        double cd    = sin(delta_0);
        double tand  = sin(delta_0);

        // Case A: theta0 = 0
        if (theta0 == 0.0) {
            if (delta_0 == 0.0 && lonpole == 90.0)
                delta_p = latpole;
            else
                delta_p = acos(sd/cp);
            if (latpole != 90.0) {
                if (fabs(latpole + delta_p) < fabs(latpole - delta_p))
                    delta_p = -delta_p;
            }
            if (lonpole == 180.0 || cd == 0.0)
                alpha_p = alpha_0;
            else
                alpha_p = alpha_0 - atan2(sp/cd, -tan(delta_p)*tand);
        }

        // Case B: theta0 != 0
        else {
            double ctheta = cos(theta0*deg2rad);
            double stheta = sin(theta0*deg2rad);
            double term1  = atan2(stheta, ctheta*cp); 
            double term2  = acos(sd/(sqrt(1.0-ctheta*ctheta*sp*sp)));
            if (term2 == 0.0)
                delta_p = term1;
            else {
                double delta_p1 = fabs((term1+term2)*rad2deg);
                double delta_p2 = fabs((term1-term2)*rad2deg);
                if (delta_p1 > 90.0 && delta_p2 > 90.0) // NO VALID SOLUTION
                    delta_p = 90.0*deg2rad;             // TODO: INVALID
                else if (delta_p1 <= 90.0 && delta_p2 > 90.0)
                    delta_p = term1 + term2;
                else if (delta_p1 > 90.0 && delta_p2 <= 90.0)
                    delta_p = term1 - term2;
                else { // Two valid solutions
                    delta_p1 = (term1+term2)*rad2deg;
                    delta_p2 = (term1-term2)*rad2deg;
                    if (fabs(latpole-delta_p1) < fabs(latpole-delta_p2))
                        delta_p = term1 + term2;
                    else
                        delta_p = term1 - term2;
                }
                if (cd == 0.0)
                    alpha_p = alpha_0;
                else {
                    double sdelt = sin(delta_p);
                    if (sdelt == 1) 
                        alpha_p = alpha_0 - phi_p - pi;
                    else {
                        if (sdelt == -1) 
                            alpha_p = alpha_0 - phi_p;
                        else {
                            alpha_p = alpha_0 - 
                                      atan2((stheta-sin(delta_p)*sd)/
                                            (cos(delta_p)*cd),
                                            sp*ctheta/cd);
                        } // endelse: sdelt != -1
                    } // endelse: sdelt != 1
                } // endelse: cd != 0
            } // endelse: term2 != 0
        } // endelse: Case B
    } // endelse: theta0 != 90

    // Set result vector
    GVector pole(alpha_p,delta_p);

    // Return pole
    return pole;
}


/***********************************************************************//**
 * @brief Get matrix for rotation between sky and native coordinates
 *
 * Requires m_npole and m_native_pole to be set.
 ***************************************************************************/
GMatrix GWcs::wcs_get_rot(void)
{
    // Allocate rotation matrix
    GMatrix r(3,3);

    // Compute useful quantities relating to reference angles
    double sp = sin(m_npole[0]*deg2rad);
    double cp = cos(m_npole[0]*deg2rad);
    double sa = sin(m_native_pole[0]);
    double ca = cos(m_native_pole[0]);
    double sd = sin(m_native_pole[1]);
    double cd = cos(m_native_pole[1]);

    // Compute rotation matrix
    r(0,0) = -sa*sp - ca*cp*sd;
    r(0,1) =  ca*sp - sa*cp*sd;
    r(0,2) =  cp*cd;
    r(1,0) =  sa*cp - ca*sp*sd;
    r(1,1) = -ca*cp - sa*sp*sd;
    r(1,2) =  sp*cd;
    r(2,0) =  ca*cd;
    r(2,1) =  sa*cd;
    r(2,2) =  sd;

    // Return rotation matrix
    return r;
}


/***********************************************************************//**
 * @brief Set standard WCS parameters
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel (first pixel is 1).
 * @param[in] crpix2 Y index of reference pixel (first pixel is 1).
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 *
 * @exception GException::wcs
 *            Unable to invert CD matrix.
 *
 * This method sets the WCS standard parameters as class members. It does
 * however not set the rotation matrices. A call to wcs_init is mandatory
 * before the projection can be used.
 *
 * @todo Implement parameter check
 ***************************************************************************/
void GWcs::wcs_set(const std::string& coords,
                   const double& crval1, const double& crval2,
                   const double& crpix1, const double& crpix2,
                   const double& cdelt1, const double& cdelt2)

{
    //TODO: Check parameters

    // Set coordinate system
    coordsys(coords);

    // Set parameters
    m_crval   = GVector(crval1, crval2);
    m_crpix   = GVector(crpix1, crpix2);
    m_cdelt   = GVector(cdelt1, cdelt2);
    m_refval  = m_crval;
    m_refpix  = m_crpix - GVector(1.0,1.0);

    // Set CD matrix corresponding to no rotation
    m_cd(0,0) = m_cdelt[0];
    m_cd(1,0) = 0.0;
    m_cd(0,1) = 0.0;
    m_cd(1,1) = m_cdelt[1];

    // Compute inverse CD matrix
    double delta = m_cd(0,0)*m_cd(1,1) - m_cd(0,1)*m_cd(1,0);
    if (delta != 0.0) {
        m_invcd(0,0) =  m_cd(1,1) / delta;
        m_invcd(0,1) = -m_cd(0,1) / delta;
        m_invcd(1,0) = -m_cd(1,0) / delta;
        m_invcd(1,1) =  m_cd(0,0) / delta;
    }
    else
        throw GException::wcs(G_WCS_SET, "Unable to invert CD matrix.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read WCS definiton from FITS header.
 *
 * @param[in] hdu FITS HDU containing the WCS definition.
 *
 * @exception GException::fits_key_not_found
 *            Unable to find required FITS header keyword.
 * @exception GException::wcs_bad_coords
 *            Coordinate system is not of valid type.
 *
 * This method reads the WCS definition from the FITS header. It calls the
 * wcs_set() method
 *
 * @todo Projection ketwords are not yet extracted from the header. Assume
 *       theta0=0.0 for the moment.
 ***************************************************************************/
void GWcs::wcs_read(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Get standard keywords
        std::string ctype1 = hdu->string("CTYPE1");
        std::string ctype2 = hdu->string("CTYPE2");
        double      crval1 = hdu->real("CRVAL1");
        double      crval2 = hdu->real("CRVAL2");
        double      crpix1 = hdu->real("CRPIX1");
        double      crpix2 = hdu->real("CRPIX2");
        double      cdelt1 = hdu->real("CDELT1");
        double      cdelt2 = hdu->real("CDELT2");

        // Determine coordinate system
        std::string coords;
        std::string xcoord = ctype1.substr(0,4);
        std::string ycoord = ctype2.substr(0,4);
        if (xcoord == "RA--" && ycoord == "DEC-") {
            m_reverse = 0;
            coords    = "EQU";
        }
        else if (xcoord == "DEC-" && ycoord == "RA--") {
            m_reverse = 1;
            coords    = "EQU";
        }
        else if (xcoord == "GLON" && ycoord == "GLAT") {
            m_reverse = 0;
            coords    = "GAL";
        }
        else if (xcoord == "GLAT" && ycoord == "GLON") {
            m_reverse = 1;
            coords    = "GAL";
        }
        else
            throw GException::wcs_bad_coords(G_WCS_READ, coordsys());

        // Set standard parameters
        wcs_set(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);

        // Get projection keywords
        //TODO

        // Initialise derived projection parameters
        wcs_init(0.0);  // theta0 = 0.0

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write WCS definiton into FITS HDU
 *
 * @param[in] hdu FITS HDU into which WCS definition will be written.
 *
 * @exception GException::wcs_bad_coords
 *            Coordinate system is not of valid type.
 ***************************************************************************/
void GWcs::wcs_write(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Set coordinate system dependent strings
        std::string c_equinox;
        std::string c_crval1;
        std::string c_crval2;
        std::string c_cdelt1;
        std::string c_cdelt2;
        std::string c_lonpole;
        std::string c_latpole;
        if (m_coordsys == 0) {
            c_equinox = "Equinox for celestial coordinate system";
            c_crval1  = "[deg] Right Ascension at reference pixel";
            c_crval2  = "[deg] Declination at reference pixel";
            c_cdelt1  = "[deg/pixel] Right Ascension increment at"
                        " reference pixel";
            c_cdelt2  = "[deg/pixel] Declination increment at"
                        " reference pixel";
            c_lonpole = "[deg] Native longitude of Celestial North Pole";
            c_latpole = "[deg] Native latitude of Celestial North Pole";
        }
        else if (m_coordsys == 1) {
            c_equinox = "Equinox for celestial coordinate system";
            c_crval1 = "[deg] Galactic Longitude at reference pixel";
            c_crval2 = "[deg] Galactic Latitude at reference pixel";
            c_cdelt1  = "[deg/pixel] Galactic Longitude increment at"
                        " reference pixel";
            c_cdelt2  = "[deg/pixel] Galactic Latitude increment at"
                        " reference pixel";
            c_lonpole = "[deg] Native longitude of Galactic North Pole";
            c_latpole = "[deg] Native latitude of Galactic North Pole";
        }
        else
            throw GException::wcs_bad_coords(G_WCS_WRITE, coordsys());

        // Set keywords
        hdu->card("EQUINOX", 2000.0,       c_equinox);
        hdu->card("CTYPE1",  wcs_crval1(), "Projection Type");
        hdu->card("CTYPE2",  wcs_crval2(), "Projection Type");
        hdu->card("CDELT1",  m_cdelt[0],   c_cdelt1);
        hdu->card("CDELT2",  m_cdelt[1],   c_cdelt2);
        hdu->card("CRPIX1",  m_crpix[0],
                  "X index of reference pixel (starting from 1)");
        hdu->card("CRPIX2",  m_crpix[1],
                  "Y index of reference pixel (starting from 1)");
        hdu->card("CRVAL1",  m_crval[0],   c_crval1);
        hdu->card("CRVAL2",  m_crval[1],   c_crval2);
        hdu->card("CROTA2",  0.0,          "[deg] Rotation Angle");
        hdu->card("LONPOLE", m_npole[0],   c_lonpole);
        hdu->card("LATPOLE", m_npole[1],   c_lonpole);
        hdu->card("PV2_1",   0.0,          "Projection parameter 1");
        hdu->card("PV2_2",   0.0,          "Projection parameter 2");

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return crval1 string
 *
 * @exception GException::wcs_invalid
 *            WCS projection not valid.
 * @exception GException::wcs_bad_coords
 *            Coordinate system is not of valid type.
 ***************************************************************************/
std::string GWcs::wcs_crval1(void) const
{
    // Initialise result
    std::string crval;

    // Check on correct type length
    if (type().length() != 3)
        throw GException::wcs_invalid(G_WCS_CRVAL1, type(),
              "3-character type required.");

    // Set coordinate system dependent value
    if (m_coordsys == 0)
        crval = "RA---" + type();
    else if (m_coordsys == 1)
        crval = "GLON-" + type();
    else
        throw GException::wcs_bad_coords(G_WCS_CRVAL1, coordsys());

    // Return
    return crval;
}


/***********************************************************************//**
 * @brief Return crval2 string
 *
 * @exception GException::wcs_invalid
 *            WCS projection not valid.
 * @exception GException::wcs_bad_coords
 *            Coordinate system is not of valid type.
 ***************************************************************************/
std::string GWcs::wcs_crval2(void) const
{
    // Initialise result
    std::string crval;

    // Check on correct type length
    if (type().length() != 3)
        throw GException::wcs_invalid(G_WCS_CRVAL2, type(),
              "3-character type required.");

    // Set coordinate system dependent value
    if (m_coordsys == 0)
        crval = "DEC--" + type();
    else if (m_coordsys == 1)
        crval = "GLAT-" + type();
    else
        throw GException::wcs_bad_coords(G_WCS_CRVAL1, coordsys());

    // Return
    return crval;
}


/***********************************************************************//**
 * @brief Print WCS information
 ***************************************************************************/
std::string GWcs::wcs_dump(void) const
{
    // Initialise result string
    std::string result;

    // Append information
    GVector native = m_native_pole*rad2deg; 
    result.append(parformat("Coordinate system")+coordsys()+"\n");
    result.append(parformat("Reference coordinate")+m_crval.print()+" deg\n");
    result.append(parformat("Reference pixel")+m_crpix.print()+"\n");
    result.append(parformat("Increment at reference")+m_cdelt.print()+" deg\n");
    result.append(parformat("Coordinate of North Pole")+m_npole.print()+" deg\n");

    // Append CD matrix
    result.append(parformat("CD matrix"));
    result.append("["+str(m_cd(0,0)));
    result.append(" "+str(m_cd(0,1)));
    result.append("]["+str(m_cd(1,0)));
    result.append(" "+str(m_cd(1,1))+"]\n");

    // Append inverse CD matrix
    result.append(parformat("Inverse CD matrix"));
    result.append("["+str(m_invcd(0,0)));
    result.append(" "+str(m_invcd(0,1)));
    result.append("]["+str(m_invcd(1,0)));
    result.append(" "+str(m_invcd(1,1))+"]\n");

    // Append coordinate of native pole
    result.append(parformat("Coordinate of native pole")+native.print());

    // Return
    return result;
}


/*==========================================================================
 =                                                                         =
 =                           GWcs private methods                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GWcs::init_members(void)
{
    // Initialise members
    m_type.clear();
    m_coordsys = 0;
    m_reverse  = 0;
    m_crval    = GVector(2);
    m_crpix    = GVector(2);
    m_cdelt    = GVector(2);
    m_npole    = GVector(180.0, 0.0);
    m_cd       = GMatrix(2,2);
    m_pv2      = GVector(21);

    // Initialise derived members
    m_theta0      = 0.0;
    m_refval      = GVector(2);
    m_refpix      = GVector(2);
    m_invcd       = GMatrix(2,2);
    m_native_pole = GVector(2);
    m_rot         = GMatrix(3,3);
    m_trot        = GMatrix(3,3);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs GWcs instance from which members should be copied
 ***************************************************************************/
void GWcs::copy_members(const GWcs& wcs)
{
    // Copy attributes
    m_type     = wcs.m_type;
    m_coordsys = wcs.m_coordsys;
    m_reverse  = wcs.m_reverse;
    m_crval    = wcs.m_crval;
    m_crpix    = wcs.m_crpix;
    m_cdelt    = wcs.m_cdelt;
    m_npole    = wcs.m_npole;
    m_cd       = wcs.m_cd;
    m_pv2      = wcs.m_pv2;

    // Copy derived attributes
    m_theta0      = wcs.m_theta0;
    m_refval      = wcs.m_refval;
    m_refpix      = wcs.m_refpix;
    m_invcd       = wcs.m_invcd;
    m_native_pole = wcs.m_native_pole;
    m_rot         = wcs.m_rot;
    m_trot        = wcs.m_trot;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcs::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Equality operator
 *
 * @param[in] a First WCS.
 * @param[in] b Second WCS.
 ***************************************************************************/
bool operator== (const GWcs &a, const GWcs &b)
{
    // Set result
    bool result = ((a.m_type     == b.m_type) &&
                   (a.m_coordsys == b.m_coordsys) &&
                   (a.m_reverse  == b.m_reverse) &&
                   (a.m_crval    == b.m_crval) &&
                   (a.m_crpix    == b.m_crpix) &&
                   (a.m_cdelt    == b.m_cdelt) &&
                   (a.m_npole    == b.m_npole) &&
                   (a.m_cd       == b.m_cd));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] a First WCS.
 * @param[in] b Second WCS.
 ***************************************************************************/
bool operator!= (const GWcs &a, const GWcs &b)
{
    // Return result
    return !(a == b);
}


/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] wcs WCS.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GWcs& wcs)
{
     // Write WCS in output stream
    os << wcs.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] wcs WCS.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GWcs& wcs)
{
    // Write WCS into logger
    log << wcs.print();

    // Return logger
    return log;
}
