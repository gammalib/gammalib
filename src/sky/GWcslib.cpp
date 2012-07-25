/***************************************************************************
 *           GWcslib.cpp  -  Virtual base class for wcslib based WCS       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @file GWcslib.cpp
 * @brief Implementation of virtual base class for wcslib based WCS
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcslib.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                                     "GWcslib::read(GFitsHDU*)"
#define G_OMEGA                                        "GWcslib::omega(int&)"
#define G_PIX2DIR                                    "GWcslib::pix2dir(int&)"
#define G_DIR2PIX                                "GWcslib::dir2pix(GSkyDir&)"
#define G_WCS_SET_CTYPE                            "GWcslib::wcs_set_ctype()"
#define G_WCS_P2S "GWcslib::wcs_s2p(int,int,double*,double*,double*,double*,"\
                                                              "double*,int*)"
#define G_WCS_S2P "GWcslib::wcs_s2p(int,int,double*,double*,double*,double*,"\
                                                              "double*,int*)"
#define G_CEL_SET                                        "GWcslib::cel_set()"
#define G_LIN_MATINV              "GWcslib::lin_matinv(std::vector<double>&,"\
                                                      "std::vector<double>&)"


/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_LIN_MATINV_FORCE_PC                             // Force PC usage

/* __ Debug definitions __________________________________________________ */
//#define G_DIR2XY_DEBUG                                      // Debug dir2xy
//#define G_XY2DIR_DEBUG                                      // Debug xy2dir

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */
const double GWcslib::UNDEFINED = 987654321.0e99;


/*==========================================================================
 =                                                                         =
 =                     GWcslib constructors/destructors                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcslib::GWcslib(void) : GWcs()
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
 * Construct standard WCS sky map from standard definition parameters.
 ***************************************************************************/
GWcslib::GWcslib(const std::string& coords,
                 const double& crval1, const double& crval2,
                 const double& crpix1, const double& crpix2,
                 const double& cdelt1, const double& cdelt2) : GWcs()
{
    // Initialise class members
    init_members();
    
    // Set standard parameters
    set_members(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct from FITS HDU table
 *
 * @param[in] hdu FITS HDU.
 ***************************************************************************/
GWcslib::GWcslib(const GFitsHDU* hdu) : GWcs()
{
    // Initialise class members
    init_members();

    // Read WCS definition from FITS HDU
    read(hdu);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param wcs World Coordinate System.
 ***************************************************************************/
GWcslib::GWcslib(const GWcslib& wcs) : GWcs(wcs)
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
GWcslib::~GWcslib(void)
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
GWcslib& GWcslib::operator= (const GWcslib& wcs)
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
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Read WCS definiton from FITS header.
 *
 * @param[in] hdu FITS HDU containing the WCS definition.
 *
 * @exception GException::wcs_bad_coords
 *            Coordinate system is not of valid type.
 *
 * This method reads the WCS definition from the FITS header.
 ***************************************************************************/
void GWcslib::read(const GFitsHDU* hdu)
{
    // Clear object
    clear();

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
            m_lng  = 0;
            m_lat  = 1;
            coords = "EQU";
        }
        else if (xcoord == "DEC-" && ycoord == "RA--") {
            m_lng  = 1;
            m_lat  = 0;
            coords = "EQU";
        }
        else if (xcoord == "GLON" && ycoord == "GLAT") {
            m_lng  = 0;
            m_lat  = 1;
            coords = "GAL";
        }
        else if (xcoord == "GLAT" && ycoord == "GLON") {
            m_lng  = 1;
            m_lat  = 0;
            coords = "GAL";
        }
        else if (xcoord == "ELON" && ycoord == "ELAT") {
            m_lng  = 0;
            m_lat  = 1;
            coords = "ECL";
        }
        else if (xcoord == "ELAT" && ycoord == "ELON") {
            m_lng  = 1;
            m_lat  = 0;
            coords = "ECL";
        }
        else if (xcoord == "HLON" && ycoord == "HLAT") {
            m_lng  = 0;
            m_lat  = 1;
            coords = "HEL";
        }
        else if (xcoord == "HLAT" && ycoord == "HLON") {
            m_lng  = 1;
            m_lat  = 0;
            coords = "HEL";
        }
        else if (xcoord == "SLON" && ycoord == "SLAT") {
            m_lng  = 0;
            m_lat  = 1;
            coords = "SGL";
        }
        else if (xcoord == "SLAT" && ycoord == "SLON") {
            m_lng  = 1;
            m_lat  = 0;
            coords = "SGL";
        }
        else {
            throw GException::wcs_bad_coords(G_READ, coordsys());
        }

        // Set standard parameters
        set_members(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write WCS definiton into FITS HDU header
 *
 * @param[in] hdu FITS HDU.
 *
 * This method writes the World Coordinate System definition into the FITS
 * HDU header.
 *
 * The method does nothing if the HDU is not valid.
 *
 * This method has been adapted from wcshdr.c::wcshdo().
 ***************************************************************************/
void GWcslib::write(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Write reference pixel coordinates
        for (int i = 0; i < m_naxis; ++i) {
            std::string keyname = "CRPIX"+str(i+1);
            std::string comment = "Pixel coordinate of reference point (starting from 1)";
            hdu->card(keyname, m_crpix.at(i), comment);
        }
        
        //TODO: Write linear transformation matrix
        
        // Write coordinate increment at reference point
        for (int i = 0; i < m_naxis; ++i) {
            std::string keyname = "CDELT"+str(i+1);
            std::string comment;
            if (m_cunit.at(i).length() > 0) {
                comment += "["+strip_whitespace(m_cunit.at(i))+"] ";
            }
            comment += "Coordinate increment at reference point";
            hdu->card(keyname, m_cdelt.at(i), comment);
        }
        
        // Write units of coordinate increment and reference value
        for (int i = 0; i < m_naxis; ++i) {
            if (m_cunit.at(i).length() > 0) {
                std::string keyname = "CUNIT"+str(i+1);
                std::string comment = "Units of coordinate increment and value";
                hdu->card(keyname, m_cunit.at(0), comment);
            }
        }
        
        // Write coordinate type
        wcs_set_ctype();
        for (int i = 0; i < m_naxis; ++i) {
            if (i == m_lng) {
                std::string keyname = "CTYPE"+str(i+1);
                hdu->card(keyname, m_ctype.at(i), m_ctype_c.at(i));
            }
            if (i == m_lat) {
                std::string keyname = "CTYPE"+str(i+1);
                hdu->card(keyname, m_ctype.at(i), m_ctype_c.at(i));
            }
            if (i == m_spec) {
                std::string keyname = "CTYPE"+str(i+1);
                hdu->card(keyname, m_ctype.at(i), m_ctype_c.at(i));
            }
        }

        // Write coordinate value at reference point
        for (int i = 0; i < m_naxis; ++i) {
            std::string keyname = "CRVAL"+str(i+1);
            std::string comment;
            if (m_cunit.at(i).length() > 0) {
                comment += "["+strip_whitespace(m_cunit.at(i))+"] ";
            }
            comment += "Coordinate value at reference point";
            hdu->card(keyname, m_crval.at(i), comment);
        }
        
        //TODO: Parameter values
        hdu->card("CROTA2", 0.0, "[deg] Rotation Angle"); // Old style, use PV instead
        //hdu->card("PV2_1",   0.0,          "Projection parameter 1");
        //hdu->card("PV2_2",   0.0,          "Projection parameter 2");
        
        // Celestial and spectral transformation parameters
        if (!undefined(m_lonpole)) {
            hdu->card("LONPOLE", m_lonpole, "[deg] Native longitude of celestial pole");
        }
        if (!undefined(m_latpole)) {
            hdu->card("LATPOLE", m_latpole, "[deg] Native latitude of celestial pole");
        }
        if (!undefined(m_restfrq)) {
            hdu->card("RESTFRQ", m_restfrq, "[Hz] Line rest frequency");
        }
        if (!undefined(m_restwav)) {
            hdu->card("RESTWAV", m_restwav, "[Hz] Line rest wavelength");
        }
        
        // Equatorial coordinate system type
        if (m_radesys.length() > 0) {
            hdu->card("RADESYS", m_radesys, "Equatorial coordinate system");
        }
        
        // Equinox of equatorial coordinate system
        if (!undefined(m_equinox)) {
            hdu->card("EQUINOX", m_equinox, "[yr] Equinox of equatorial coordinates");
        }
        
        //TODO: Reference frame of spectral coordinates
        
        //TODO: Reference frame of spectral observation
        
        //TODO: Observer's velocity towards source
        
        //TODO: Reference frame of source redshift
        
        //TODO: Redshift of the source
        
        //TODO: Observatory coordinates
        
        //TODO: MJD of observation
        
        //TODO: MJD mid-observation time
        
        //TODO: ISO-8601 date corresponding to MJD-OBS
        
        //TODO: ISO-8601 date corresponding to MJD-AVG

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel in units of steradians
 *
 * @param[in] pix Pixel index.
 *
 * This dummy method throws an error when called as pixel indices are not
 * implemented for a HealPix grid. Maybe we should drop this method and
 * implement a general GSkyPixel class tha also handles HealPix?
 *
 * @todo Think about implementation of sky pixel for HealPix data. We may
 *       use only the first argument, and even carry a usage flag in
 *       GSkyPixel
 ***************************************************************************/
double GWcslib::omega(const int& pix) const
{
    // Set error message
    std::string message = "Method not defined for general WCS projection.";

    // Throw error
    throw GException::wcs(G_OMEGA, message);

    // Return
    return 0.0;
}


/***********************************************************************//**
 * @brief Returns solid angle of pixel in units of steradians
 *
 * @param[in] pix Pixel index (x,y)
 *
 * Estimate solid angles of pixel by compuing the coordinates in the 4 pixel
 * corners. The surface is computed using a cartesian approximation:
 *           a
 *     1-----5-----2                 a+b
 *    /      | h    \    where A = h ---
 *   4-------6-------3                2
 *           b
 * This is a brute force technique that works sufficiently well for non-
 * rotated sky maps. Something more intelligent should be implemented in
 * the future.
 ***************************************************************************/
double GWcslib::omega(const GSkyPixel& pix) const
{
    // Get the sky directions of the 6 points
    GSkyDir dir1 = xy2dir(GSkyPixel(pix.x()-0.5, pix.y()-0.5));
    GSkyDir dir2 = xy2dir(GSkyPixel(pix.x()+0.5, pix.y()-0.5));
    GSkyDir dir3 = xy2dir(GSkyPixel(pix.x()+0.5, pix.y()+0.5));
    GSkyDir dir4 = xy2dir(GSkyPixel(pix.x()-0.5, pix.y()+0.5));
    GSkyDir dir5 = xy2dir(GSkyPixel(pix.x(), pix.y()-0.5));
    GSkyDir dir6 = xy2dir(GSkyPixel(pix.x(), pix.y()+0.5));

    // Compute distances between sky directions
    double a = dir1.dist(dir2);
    double b = dir3.dist(dir4);
    double h = dir5.dist(dir6);

    // Compute solid angle
    double omega = 0.5*(h*(a+b));

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
GSkyDir GWcslib::pix2dir(const int& pix) const
{
    // Set error message
    std::string message = "Method not defined for general WCS projection.";

    // Throw error
    throw GException::wcs(G_PIX2DIR, message);

    // Define direction
    GSkyDir dir;

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
int GWcslib::dir2pix(const GSkyDir& dir) const
{
    // Set error message
    std::string message = "Method not defined for general WCS projection.";

    // Throw error
    throw GException::wcs(G_DIR2PIX, message);

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
GSkyDir GWcslib::xy2dir(const GSkyPixel& pix) const
{
    // Allocate memory for transformation
    double pixcrd[2];
    double imgcrd[2];
    double phi;
    double theta;
    double world[2];
    int    stat;
    
    // Set sky pixel. We have to add 1.0 here as the WCS pixel reference
    // (CRPIX) starts from one while GSkyPixel starts from 0.
    pixcrd[0] = pix.x() + 1.0;
    pixcrd[1] = pix.y() + 1.0;
    
    // Trasform pixel-to-world coordinate
    wcs_p2s(1, 2, pixcrd, imgcrd, &phi, &theta, world, &stat);

    // Set sky direction
    GSkyDir dir;
    if (m_coordsys == 0) {
        dir.radec_deg(world[0], world[1]);
    }
    else {
        dir.lb_deg(world[0], world[1]);
    }

    // Debug: Dump transformation steps
    #if defined(G_XY2DIR_DEBUG)
    std::cout << "xy2dir: pixel=" << pix
              << " (x,y)=(" << pixcrd[0] << "," << pixcrd[1] << ")"
              << " (phi,theta)=(" << phi << "," << theta << ")"
              << " (lng,lat)=(" << world[0] << "," << world[1] << ")" << std::endl;
    #endif

    // Return
    return dir;
}


/***********************************************************************//**
 * @brief Returns pixel of sky direction
 *
 * @param[in] dir Sky direction.
 *
 * Note that GSkyPixel indices start from 0 while the WCS pixel reference
 * starts from 1.
 ***************************************************************************/
GSkyPixel GWcslib::dir2xy(const GSkyDir& dir) const
{
    // Allocate memory for transformation
    double pixcrd[2];
    double imgcrd[2];
    double phi;
    double theta;
    double world[2];
    int    stat;
    
    // Set world coordinate
    if (m_coordsys == 0) {
        world[0] = dir.ra_deg();
        world[1] = dir.dec_deg();
    }
    else {
        world[0] = dir.l_deg();
        world[1] = dir.b_deg();
    }
    
    // Transform world-to-pixel coordinate
    wcs_s2p(1, 2, world, &phi, &theta, imgcrd, pixcrd, &stat);
    
    // Set sky pixel. We have to subtract 1 here as GSkyPixel starts from
    // zero while the WCS reference (CRPIX) starts from one.
    GSkyPixel pixel(pixcrd[0]-1.0, pixcrd[1]-1.0);

    // Debug: Dump transformation steps
    #if defined(G_DIR2XY_DEBUG)
    std::cout << "dir2xy: dir=" << dir
              << " (lng,lat)=(" << world[0] << "," << world[1] << ")"
              << " (phi,theta)=(" << phi << "," << theta << ")"
              << " (x,y)=" << pixel << std::endl;
    #endif

    // Return sky pixel
    return pixel;
}


/***********************************************************************//**
 * @brief Set World Coordinate System parameters
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel (first pixel is 1).
 * @param[in] crpix2 Y index of reference pixel (first pixel is 1).
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 *
 * This method sets the WCS parameters.
 ***************************************************************************/
void GWcslib::set(const std::string& coords,
                  const double& crval1, const double& crval2,
                  const double& crpix1, const double& crpix2,
                  const double& cdelt1, const double& cdelt2)

{
    // Clear any existing information
    clear();
    
    // Set standard parameters
    set_members(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2);
    
    // Setup WCS derived parameters
    wcs_set();
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 *
 * Code adapted from wcs.c::wcsini().
 ***************************************************************************/
void GWcslib::init_members(void)
{
    // Initialise World Coordinate System with 0 axes
    wcs_ini(0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs GWcs instance from which members should be copied
 ***************************************************************************/
void GWcslib::copy_members(const GWcslib& wcs)
{
    // Copy WCS attributes
    m_wcsset  = wcs.m_wcsset;
    m_naxis   = wcs.m_naxis;
    m_crval   = wcs.m_crval;
    m_cunit   = wcs.m_cunit;
    m_ctype   = wcs.m_ctype;
    m_ctype_c = wcs.m_ctype_c;
    m_lonpole = wcs.m_lonpole;
    m_latpole = wcs.m_latpole;
    m_restfrq = wcs.m_restfrq;
    m_restwav = wcs.m_restwav;
    m_radesys = wcs.m_radesys;
    m_equinox = wcs.m_equinox;
    m_cd      = wcs.m_cd;
    m_crota   = wcs.m_crota;
    m_lng     = wcs.m_lng;
    m_lat     = wcs.m_lat;
    m_spec    = wcs.m_spec;

    // Copy linear transformation parameters
    m_linset  = wcs.m_linset;
    m_unity   = wcs.m_unity;
    m_crpix   = wcs.m_crpix;
    m_pc      = wcs.m_pc;
    m_cdelt   = wcs.m_cdelt;
    m_piximg  = wcs.m_piximg;
    m_imgpix  = wcs.m_imgpix;

    // Copy celestial transformation parameters
    m_celset   = wcs.m_celset;
    m_offset   = wcs.m_offset;
    m_phi0     = wcs.m_phi0;
    m_theta0   = wcs.m_theta0;
    m_ref[0]   = wcs.m_ref[0];
    m_ref[1]   = wcs.m_ref[1];
    m_ref[2]   = wcs.m_ref[2];
    m_ref[3]   = wcs.m_ref[3];
    m_euler[0] = wcs.m_euler[0];
    m_euler[1] = wcs.m_euler[1];
    m_euler[2] = wcs.m_euler[2];
    m_euler[3] = wcs.m_euler[3];
    m_euler[4] = wcs.m_euler[4];
    m_latpreq  = wcs.m_latpreq;
    m_isolat   = wcs.m_isolat;

    // Copy projection parameters
    m_prjset = wcs.m_prjset;
    m_r0     = wcs.m_r0;
    m_bounds = wcs.m_bounds;
    m_x0     = wcs.m_x0;
    m_y0     = wcs.m_y0;
    m_w      = wcs.m_w;
    for (int i = 0; i < PVN; ++i)
        m_pv[i] = wcs.m_pv[i];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcslib::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set World Coordinate System parameters
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel (first pixel is 1).
 * @param[in] crpix2 Y index of reference pixel (first pixel is 1).
 * @param[in] cdelt1 Increment in x direction at reference pixel (deg).
 * @param[in] cdelt2 Increment in y direction at reference pixel (deg).
 *
 * This method sets the WCS parameters. It does not call wcs_set(), however,
 * as wcs_set() may depend on the availability of derived class methods which
 * can not be used in a base class constructor.
 *
 * @todo Implement parameter validity check
 ***************************************************************************/
void GWcslib::set_members(const std::string& coords,
                          const double& crval1, const double& crval2,
                          const double& crpix1, const double& crpix2,
                          const double& cdelt1, const double& cdelt2)

{
    //TODO: Check parameters

    // Initialise WCS
    wcs_ini(2);

    // Set coordinate system
    coordsys(coords);

    // Set World Coordinate parameters
    m_lng      = 0;
    m_lat      = 1;
    m_crpix[0] = crpix1;
    m_crpix[1] = crpix2;
    m_cdelt[0] = cdelt1;
    m_cdelt[1] = cdelt2;
    m_crval[0] = crval1;
    m_crval[1] = crval2;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns true if argument is identical
 *
 * @param[in] wcs Pointer to World Coordinate System
 *
 * This method is a helper for the World Coordinate Comparison friends. It
 * does not compare derived quantities as those may not have been initialised
 * in both objects.
 ***************************************************************************/
bool GWcslib::compare(const GWcs& wcs) const
{
    // Initialise result
    bool result = false;
    
    // Continue only we compare to a GWcslib object
    const GWcslib* ptr = dynamic_cast<const GWcslib*>(&wcs);
    if (ptr != NULL) {
    
        // Perform comparion of standard (non derived) parameters
        result = ((code()     == ptr->code())     &&
                  (m_coordsys == ptr->m_coordsys) &&
                  (m_lng      == ptr->m_lng)      &&
                  (m_lat      == ptr->m_lat)      &&
                  (m_crval    == ptr->m_crval)    &&
                  (m_crpix    == ptr->m_crpix)    &&
                  (m_cdelt    == ptr->m_cdelt));
        
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =              World Coordinate Methods adapted from wcslib               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise World Coordinate System
 *
 * @param[in] naxis Number of axes
 *
 * This method initialises the World Coordinate System information. The
 * method has been adapted from the wcslib function wcs.c::wcsini. In
 * contrast to the wcslib function, however, this method accepts naxis=0.
 * In this case, all vectors are cleared.
 ***************************************************************************/
void GWcslib::wcs_ini(int naxis)
{
    // Signal that WCS information has not been set so far. This will be done
    // in wcs_set() upon request
    m_wcsset = false;
    
    // Clear vectors
    m_crpix.clear();
    m_pc.clear();
    m_cdelt.clear();
    m_crval.clear();
    m_cunit.clear();
    m_ctype.clear();
    m_ctype_c.clear();
    m_cd.clear();
    m_crota.clear();

    // Set number of axes
    m_naxis = naxis;
    
    // Set defaults for the linear transformation (sets m_Crpix, m_pc, m_Cdelt)
    lin_ini(naxis);
    
    // Set default axes parameters
    if (naxis > 0) {
        for (int i = 0; i < naxis; ++i) {
            m_crval.push_back(0.0);
            m_cunit.push_back("");
            m_ctype.push_back("");
            m_ctype_c.push_back("");
        }
    }
    
    // Set defaults for the celestial transformation parameters
    m_lonpole = UNDEFINED;
    m_latpole = +90.0;
    
    // Set defaults for the spectral transformation parameters
    m_restfrq = UNDEFINED; // 0.0 in wcslib
    m_restwav = UNDEFINED; // 0.0 in wcslib
    
    //TODO: Default parameter values
    
    // Defaults for alternate linear transformations
    if (naxis > 0) {
        for (int i = 0; i < naxis; ++i) {
            for (int j = 0; j < naxis; ++j) {
                m_cd.push_back(0.0);
            }
            m_crota.push_back(0.0);
        }
    }
    
    //TODO: Defaults for auxiliary coordinate system information
    m_radesys.clear();
    m_equinox = UNDEFINED;
    
    // Reset derived values
    m_lng  = -1;
    m_lat  = -1;
    m_spec = -1;
    
    // Initialise celestial transformation parameters
    cel_ini();

    // Initialise spectral transformation parameters
    spc_ini();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup of World Coordinate System
 *
 * This method sets up the World Coordinate System information. In particular
 * it:
 * - initialises the celestial parameters (call to cel_ini)
 * - sets the m_ref information (reference value and native pole)
 * - sets the celestial parameters (call to cel_set)
 * - sets the linear transformation (call to lin_set)
 * The method has been adapted from the wcslib function wcs.c::wcsset.
 *
 * @todo Determine axis types from CTYPEia
 * @todo Convert to canonical units
 * @todo Do we have PVi_ma keyvalues?
 * @todo Do simple alias translations
 * @todo Update PVi_ma keyvalues
 * @todo Non-linear spectral axis present?
 * @todo Tabular axes present?
 ***************************************************************************/
void GWcslib::wcs_set(void) const
{
    //TODO: Determine axis types from CTYPEia
    
    //TODO: Convert to canonical units
    
    // Non-linear celestial axes present?
    if (m_lng >= 0) {
        
        // Initialise celestial parameters
        cel_ini();
        
        // Set CRVALia, LONPOLEa, and LATPOLEa keyvalues
        m_ref[0] = m_crval[m_lng]; // Longitude reference value
        m_ref[1] = m_crval[m_lat]; // Latitude reference value
        m_ref[2] = m_lonpole;      // LONPOLE
        m_ref[3] = m_latpole;      // LATPOLE
        
        //TODO: Do we have PVi_ma keyvalues?
        
        //TODO: Do simple alias translations
        
        // Initialize the celestial transformation routines
        m_r0 = 0.0; // Forces initialisation
        cel_set();
        
        // Update LONPOLE and LATPOLE 
        m_lonpole = m_ref[2];
        m_latpole = m_ref[3];
        
        //TODO: Update PVi_ma keyvalues
        
    } // endif: non-linear celestial axes were present
    
    //TODO: Non-linear spectral axis present?
    
    //TODO: Tabular axes present?
    
    // Initialize the linear transformation
    lin_set();
    
    // Signal that WCS is set
    m_wcsset = true;
    
    // Return
    return;
}    


/***********************************************************************//**
 * @brief Set CTYPEa keywords
 *
 * @exception GException::wcs_invalid
 *            WCS projection not valid.
 * @exception GException::wcs_bad_coords
 *            Coordinate system is not of valid type.
 *
 * This method has been inspired by code from wcshdr.c::wcshdo.
 ***************************************************************************/
void GWcslib::wcs_set_ctype(void) const
{
    // Check on correct type length
    if (code().length() != 3)
        throw GException::wcs_invalid(G_WCS_SET_CTYPE, code(),
              "3-character type required.");
    
    // Set longitude keyword
    if (m_lng >= 0) {
    
        // Set coordinate system
        if (m_coordsys == 0) {
            m_ctype.at(m_lng) = "RA---" + code();
            m_ctype_c.at(m_lng) = "Right ascension, ";
        }
        else if (m_coordsys == 1) {
            m_ctype.at(m_lng) = "GLON-" + code();
            m_ctype_c.at(m_lng) = "Galactic longitude, ";
        }
        else if (m_coordsys == 2) {
            m_ctype.at(m_lng) = "ELON-" + code();
            m_ctype_c.at(m_lng) = "Ecliptic longitude, ";
        }
        else if (m_coordsys == 3) {
            m_ctype.at(m_lng) = "HLON-" + code();
            m_ctype_c.at(m_lng) = "Helioecliptic longitude, ";
        }
        else if (m_coordsys == 4) {
            m_ctype.at(m_lng) = "SLON-" + code();
            m_ctype_c.at(m_lng) = "Supergalactic longitude, ";
        }
        else {
            throw GException::wcs_bad_coords(G_WCS_SET_CTYPE, coordsys());
        }
        
        // Add projection name to comment
        m_ctype_c.at(m_lng).append(name());
        m_ctype_c.at(m_lng).append(" projection");
    }

    // Set latitude keyword
    if (m_lat >= 0) {
    
        // Set coordinate system
        if (m_coordsys == 0) {
            m_ctype.at(m_lat) = "DEC--" + code();
            m_ctype_c.at(m_lat) = "Declination, ";
        }
        else if (m_coordsys == 1) {
            m_ctype.at(m_lat) = "GLAT-" + code();
            m_ctype_c.at(m_lat) = "Galactic latitude, ";
        }
        else if (m_coordsys == 2) {
            m_ctype.at(m_lat) = "ELAT-" + code();
            m_ctype_c.at(m_lat) = "Ecliptic latitude, ";
        }
        else if (m_coordsys == 3) {
            m_ctype.at(m_lat) = "HLAT-" + code();
            m_ctype_c.at(m_lat) = "Helioecliptic latitude, ";
        }
        else if (m_coordsys == 4) {
            m_ctype.at(m_lat) = "SLAT-" + code();
            m_ctype_c.at(m_lat) = "Supergalactic latitude, ";
        }
        else {
            throw GException::wcs_bad_coords(G_WCS_SET_CTYPE, coordsys());
        }

        // Add projection name to comment
        m_ctype_c.at(m_lng).append(name());
        m_ctype_c.at(m_lng).append(" projection");
    }

    //TODO: Set spectral keyword
    if (m_spec >= 0) {
        m_ctype.at(m_spec) = "NONE";
        m_ctype_c.at(m_spec) = "Not yet implemented";
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Pixel-to-world transformation
 *
 * @param[in] ncoord Number of coordinates.
 * @param[in] nelem Vector length of each coordinate (>=m_naxis).
 * @param[in] pixcrd Array [ncoord][nelem] of pixel coordinates.
 * @param[out] imgcrd Array [ncoord][nelem] of intermediate world coordinates.
 * @param[out] phi Array [ncoord] of longitudes in native coordinate system.
 * @param[out] theta Array [ncoord] of latitudes in native coordinate system.
 * @param[out] world Array [ncoord][nelem] of world coordinates.
 * @param[out] stat Array [ncoord] of pixel coordinates validity.
 *
 * @exception GException::wcs_invalid_parameter
 *            Invalid input parameters provided
 *
 * This method transforms pixel coordinates to world coordinates. The method
 * has been adapted from wcs.c::wcsp2s(). Note that this method is extremely
 * simplified with respect to wcslib, but it does the job for now.
 *
 * For celestial axes, imgcrd[][wcs.lng] and imgcrd[][wcs.lat] are the
 * projected x-, and y-coordinates in pseudo "degrees" and world[][wcs.lng]
 * and world[][wcs.lat] are the celestial longitude and latitude [deg]
 * For spectral axes, imgcrd[][wcs.spec] is the intermediate spectral
 * coordinate, in SI units and world[][wcs.spec] is ...
 *
 * @todo Check for constant x and/or y to speed-up computations
 * @todo Zero the unused world coordinate elements
 ***************************************************************************/
void GWcslib::wcs_p2s(int ncoord, int nelem, const double* pixcrd, double* imgcrd,
                      double* phi, double* theta, double* world, int* stat) const
{
    // Initialize if required
    if (!m_wcsset)
        wcs_set();
    
    // Sanity check
    if (ncoord < 1 || (ncoord > 1 && nelem < m_naxis)) {
        std::string message;
        if (ncoord < 1) {
            message = "ncoord="+str(ncoord)+", >0 required.";
        }
        else {
            message = "nelem="+str(nelem)+", >="+str(m_naxis)+" required.";
        }
        throw GException::wcs_invalid_parameter(G_WCS_P2S, message);
    }
    
    // Apply pixel-to-world linear transformation
    lin_p2x(ncoord, nelem, pixcrd, imgcrd);
    
    //TODO: Check for constant x and/or y
    int nx = ncoord;
    int ny = 0;
        
    // Transform projection plane coordinates to celestial coordinates
    cel_x2s(nx, ny, nelem, nelem, imgcrd+m_lng, imgcrd+m_lat, phi, theta,
            world+m_lng, world+m_lat, stat);

    //TODO: Zero the unused world coordinate elements
    
    // Return
    return;
}    


/***********************************************************************//**
 * @brief World-to-pixel transformation
 *
 * @param[in] ncoord Number of coordinates.
 * @param[in] nelem Vector length of each coordinate (>=m_naxis).
 * @param[in] world Array [ncoord][nelem] of world coordinates.
 * @param[out] phi Array [ncoord] of longitudes in native coordinate system.
 * @param[out] theta Array [ncoord] of latitudes in native coordinate system.
 * @param[out] imgcrd Array [ncoord][nelem] of intermediate world coordinates.
 * @param[out] pixcrd Array [ncoord][nelem] of pixel coordinates.
 * @param[out] stat Array [ncoord] of pixel coordinates validity.
 *
 * @exception GException::wcs_invalid_parameter
 *            Invalid input parameters provided
 *
 * This method transforms world coordinates to pixel coordinates. The method
 * has been adapted from wcs.c::wcss2p(). Note that this method is extremely
 * simplified with respect to wcslib, but it does the job for now.
 *
 * For celestial axes, imgcrd[][wcs.lng] and imgcrd[][wcs.lat] are the
 * projected x-, and y-coordinates in pseudo "degrees" and world[][wcs.lng]
 * and world[][wcs.lat] are the celestial longitude and latitude [deg]
 * For spectral axes, imgcrd[][wcs.spec] is the intermediate spectral
 * coordinate, in SI units and world[][wcs.spec] is ...
 *
 * @todo Check for constant x and/or y to speed-up computations
 * @todo Zero the unused world coordinate elements
 ***************************************************************************/
void GWcslib::wcs_s2p(int ncoord, int nelem, const double* world,
                      double* phi, double* theta,  double* imgcrd,
                      double* pixcrd, int* stat) const
{
    // Initialize if required
    if (!m_wcsset)
        wcs_set();
    
    // Sanity check
    if (ncoord < 1 || (ncoord > 1 && nelem < m_naxis)) {
        std::string message;
        if (ncoord < 1) {
            message = "ncoord="+str(ncoord)+", >0 required.";
        }
        else {
            message = "nelem="+str(nelem)+", >="+str(m_naxis)+" required.";
        }
        throw GException::wcs_invalid_parameter(G_WCS_S2P, message);
    }

    //TODO: Check for constant x and/or y
    int nlng = ncoord;
    int nlat = 0;

    // Transform celestial coordinates to projection plane coordinates
    cel_s2x(nlng, nlat, nelem, nelem, world+m_lng, world+m_lat,
            phi, theta, imgcrd+m_lng, imgcrd+m_lat, stat);

    //TODO: Zero the unused world coordinate elements

    // Apply world-to-pixel linear transformation
    lin_x2p(ncoord, nelem, imgcrd, pixcrd);
    
    // Return
    return;
}    


/***********************************************************************//**
 * @brief Print WCS information
 ***************************************************************************/
std::string GWcslib::wcs_print(void) const
{
    // Initialise result string
    std::string result;

    // Append World Coordinate parameters
    result.append(parformat("Number of axes")+str(m_naxis)+"\n");
    result.append(parformat("Longitude axis")+str(m_lng)+"\n");
    result.append(parformat("Latitude axis")+str(m_lat)+"\n");
    result.append(parformat("Spectral axis")+str(m_spec)+"\n");
    
    // Append coordinates
    result.append(parformat("Reference coordinate")+"(");
    for (int i = 0; i < m_crval.size(); ++i) {
        if (i > 0) {
            result.append(", ");
        }
        result.append(str(m_crval[i]));
        if (m_cunit[i].length() > 0) {
            result.append(" "+m_cunit[i]);
        }
            
    }
    result.append(")\n");
    result.append(parformat("Reference pixel")+"(");
    for (int i = 0; i < m_crpix.size(); ++i) {
        if (i > 0) {
            result.append(", ");
        }
        result.append(str(m_crpix[i]));
    }
    result.append(")\n");
    result.append(parformat("Increment at reference")+"(");
    for (int i = 0; i < m_cdelt.size(); ++i) {
        if (i > 0) {
            result.append(", ");
        }
        result.append(str(m_cdelt[i]));
        if (m_cunit[i].length() > 0) {
            result.append(" "+m_cunit[i]);
        }
    }
    result.append(")\n");
        
    // Append origin
    result.append(parformat("(Phi_0, Theta_0)")+"(");
    result.append(wcs_print_value(m_phi0)+", ");
    result.append(wcs_print_value(m_theta0)+") deg\n");
    
    // Append native pole
    result.append(parformat("(Phi_p, Theta_p)")+"(");
    result.append(wcs_print_value(m_lonpole)+", ");
    result.append(wcs_print_value(m_latpole)+") deg\n");

    // Append LATPOLEa keyword usage
    result.append(parformat("LATPOLE keyword usage"));
    switch (m_latpreq) {
    case 0:
        result.append("Not used. Theta_p determined uniquely by"
                      " CRVALia and LONPOLEa keywords.\n");
        break;
    case 1:
        result.append("Required to select between two valid solutions"
                      " of Theta_p.\n");
        break;
    case 2:
        result.append("Theta_p was specified solely by LATPOLE.\n");
        break;
    default:
        result.append("UNDEFINED\n");
        break;
    }

    // Append celestial transformation parameters
    result.append(parformat("Reference vector (m_ref)")+"(");
    for (int k = 0; k < 4; ++k) {
        if (k > 0) {
            result.append(", ");
        }
        result.append(wcs_print_value(m_ref[k]));
    }
    result.append(") deg\n");

    // Append Euler angles
    result.append(parformat("Euler angles")+"(");
    for (int k = 0; k < 5; ++k) {
        if (k > 0) {
            result.append(", ");
        }
        result.append(str(m_euler[k]));
        if (k < 3) {
            result.append(" deg");
        }
    }
    result.append(")\n");
    
    // Append latitude preservement flag
    if (m_isolat) {
        result.append(parformat("Latitude preserved")+"True\n");
    }
    else {
        result.append(parformat("Latitude preserved")+"False\n");
    }

    // Append linear transformation parameters
    if (m_unity) {
        result.append(parformat("Unity PC matrix")+"True\n");
    }
    else {
        result.append(parformat("Unity PC matrix")+"False\n");
    }
    result.append(parformat("Pixel-to-image trafo")+"(");
    for (int k = 0; k < m_piximg.size(); ++k) {
        if (k > 0) {
            result.append(", ");
        }
        result.append(str(m_piximg[k]));
    }
    result.append(")\n");
    result.append(parformat("Image-to-pixel trafo")+"(");
    for (int k = 0; k < m_imgpix.size(); ++k) {
        if (k > 0) {
            result.append(", ");
        }
        result.append(str(m_imgpix[k]));
    }
    result.append(")\n");
    
    // Append coordinate system
    result.append(parformat("Coodinate system")+coordsys()+"\n");
    
    // Append projection parameters
    result.append(parformat("Projection code")+code()+"\n");
    result.append(parformat("Projection name")+name()+"\n");
    result.append(parformat("Radius of the gen. sphere")+str(m_r0)+" deg\n");

    // Append boundary checking information
    if (m_bounds) {
        result.append(parformat("Strict bounds checking")+"True\n");
    }
    else {
        result.append(parformat("Strict bounds checking")+"False\n");
    }

    // Append fiducial offset information
    if (m_offset) {
        result.append(parformat("Use fiducial offset")+"True\n");
    }
    else {
        result.append(parformat("Use fiducial offset")+"False\n");
    }
    result.append(parformat("Fiducial offset")+"("+str(m_x0)+", "+str(m_y0)+")\n");
    
    // Append spectral transformation parameters
    result.append(parformat("Rest frequency")+wcs_print_value(m_restfrq)+"\n");
    result.append(parformat("Rest wavelength")+wcs_print_value(m_restwav));

    // Return
    return result;
}


/***********************************************************************//**
 * @brief Helper function for value printing
 *
 * @param[in] value Double precision value
 *
 * This helper function either prints a double precision value or the
 * word 'UNDEFINED', depending on whether the value is defined or not.
 ***************************************************************************/
std::string GWcslib::wcs_print_value(const double& value) const
{
    // Initialise result
    std::string result;
    
    // Depend value dependent string
    if (undefined(value)) {
        result.append("UNDEFINED");
    }
    else {
        result.append(str(value));
    }
    
    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =         Celestial transformation methods adapted from wcslib            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise celestial transformation parameters
 *
 * Code adapted from cel.c::celini().
 ***************************************************************************/
void GWcslib::cel_ini(void) const
{
    // Initialise parameters
    m_celset  = false;
    m_offset  = false;
    m_phi0    = UNDEFINED;
    m_theta0  = UNDEFINED;
    m_ref[0]  = 0.0;
    m_ref[1]  = 0.0;
    m_ref[2]  = UNDEFINED;
    m_ref[3]  = +90.0;
    m_latpreq = -1;
    m_isolat  = false;
    
    // Clear Euler angles
    for (int k = 0; k < 5; ++k) {
        m_euler[k] = 0.0;
    }
    
    // Initialise projection parameters
    prj_ini();
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup of celestial transformation
 *
 * @exception GException::wcs_invalid_parameter
 *            Invalid World Coordinate System parameter encountered.
 *
 * This method sets up the celestial transformation information. The method
 * has been adapted from the wcslib function cel.c::celset. It:
 * - sets the projection (using prj_set)
 * - computes the native longitude and latitude and stores the results 
 *   in m_ref[2] and m_ref[3]
 * - computes the Euler angles for the celestial transformation and stores
 *   the result in m_euler
 * - it sets the m_latpreq (0=, 1=two solutions)
 * - it sets the m_isolat flag
 *
 * The m_latpreq member is for informational purposes and indicates how the
 * LATPOLEa keyword was used
 * - 0: Not required, theta_p (== delta_p) was determined uniquely by the
 *      CRVALia and LONPOLEa keywords.
 * - 1: Required to select between two valid solutions of theta_p.
 * - 2: theta_p was specified solely by LATPOLEa.
 *
 * The m_isolat flag is true if the spherical rotation preserves the
 * magnitude of the latitude, which occurs if the axes of the native and
 * celestial coordinates are coincident. It signals an opportunity to cache
 * intermediate calculations common to all elements in a vector computation.
 ***************************************************************************/
void GWcslib::cel_set(void) const
{
    // Set tolerance
    const double tol = 1.0e-10;
    
    // If no fiducial offsets are requested then make sure that (phi0,theta0)
    // are undefined
    if (!m_offset) {
        m_phi0   = UNDEFINED;
        m_theta0 = UNDEFINED;
    }
    
    // Setup projection. This method will set (phi0,theta0) if they were
    // undefined.
    prj_set();
    
    // Initialise celestial parameters
    double lng0 = m_ref[0];
    double lat0 = m_ref[1];
    double phip = m_ref[2];
    double latp = m_ref[3];
    double lngp = 0.0;

    // Set default for native longitude of the celestial pole?
    if (undefined(phip) || phip == 999.0) {
        phip  = (lat0 < m_theta0) ? 180.0 : 0.0;
        phip += m_phi0;
        if (phip < -180.0) {
            phip += 360.0;
        }
        else if (phip > 180.0) {
            phip -= 360.0;
        }
        m_ref[2] = phip;
    }

    // Initialise
    m_latpreq = 0;
    
    // Compute celestial coordinates of the native pole
    // Fiducial point at the native pole?
    if (m_theta0 == 90.0) {
        lngp = lng0;
        latp = lat0;
    }
    
    // ... otherwise fiducial point is away from the native pole
    else {
    
        // Compute sine and cosine of lat0 and theta0
        double slat0;
        double clat0;
        double sthe0;
        double cthe0;
        sincosd(lat0,     &slat0, &clat0);
        sincosd(m_theta0, &sthe0, &cthe0);

        // ...
        double sphip;
        double cphip;
        double u;
        double v;
        if (phip == m_phi0) {
            sphip = 0.0;
            cphip = 1.0;
            u     = m_theta0;
            v     = 90.0 - lat0;
        }
        
        // ...
        else {
            
            // Compute sine and cosine of Phi_p - Phi_0
            sincosd(phip - m_phi0, &sphip, &cphip);

            // ...
            double x = cthe0 * cphip;
            double y = sthe0;
            double z = std::sqrt(x*x + y*y);
            if (z == 0.0) {
                
                // Check of an invalid coordinate transformation parameter
                // has been encountered. Since z=0 we exlect sin(lat0)=0.
                if (slat0 != 0.0) {
                    std::string message = "sin(lat0)="+str(slat0)+", expected 0.";
                    throw GException::wcs_invalid_parameter(G_CEL_SET, message);
                }

                // latp determined solely by LATPOLEa in this case
                m_latpreq = 2;
                if (latp > 90.0) {
                    latp = 90.0;
                }
                else if (latp < -90.0) {
                    latp = -90.0;
                }
                
            }
            
            // ... otherwise ...
            else {
                double slz = slat0/z;
                if (std::abs(slz) > 1.0) {
                    if ((std::abs(slz) - 1.0) < tol) {
                        if (slz > 0.0) {
                            slz = 1.0;
                        }
                        else {
                            slz = -1.0;
                        }
                    }
                    else {
                        std::string message;
                        message  = "abs(slz)-1 >= "+str(std::abs(slz) - 1.0);
                        message += +", expected <"+str(tol)+".";
                        throw GException::wcs_invalid_parameter(G_CEL_SET, message);
                    }
                }
                u = atan2d(y,x);
                v = acosd(slz);
            } // endelse: z != 0
        } // endelse: ...

        // ...
        if (m_latpreq == 0) {
            double latp1 = u + v;
            if (latp1 > 180.0) {
                latp1 -= 360.0;
            }
            else if (latp1 < -180.0) {
                latp1 += 360.0;
            }

            // ...
            double latp2 = u - v;
            if (latp2 > 180.0) {
                latp2 -= 360.0;
            }
            else if (latp2 < -180.0) {
                latp2 += 360.0;
            }

            // Are there two valid solutions for latp?
            if (std::abs(latp1) < 90.0+tol && std::abs(latp2) < 90.0+tol) {
                m_latpreq = 1;
            }

            // Select solution
            if (std::abs(latp-latp1) < std::abs(latp-latp2)) {
                latp = (std::abs(latp1) < 90.0+tol) ? latp1 : latp2;
            }
            else {
                latp = (std::abs(latp2) < 90.0+tol) ? latp2 : latp1;
            }

            // Account for rounding error
            if (std::abs(latp) < 90.0+tol) {
                if (latp > 90.0) {
                    latp =  90.0;
                }
                else if (latp < -90.0) {
                    latp = -90.0;
                }
            }
        } // endif: ...

        // ...
        double z = cosd(latp) * clat0;
        if (std::abs(z) < tol) {
        
            // Celestial pole at the fiducial point
            if (std::abs(clat0) < tol) {
                lngp = lng0;
            }

            // Celestial north pole at the native pole
            else if (latp > 0.0) {
                lngp = lng0 + phip - m_phi0 - 180.0;
            }

            // Celestial south pole at the native pole
            else {
                lngp = lng0 - phip + m_phi0;
            }

        }
        
        // ...
        else {
            double x = (sthe0 - sind(latp)*slat0)/z;
            double y =  sphip*cthe0/clat0;
            if (x == 0.0 && y == 0.0) {
                std::string message;
                message  = "x=0 and y=0";
                message += +", expected at least that one differs from 0.";
                throw GException::wcs_invalid_parameter(G_CEL_SET, message);
            }
            lngp = lng0 - atan2d(y,x);
        }

        // Make celestial longitude at the pole the same sign as at the
        // fiducial point
        if (lng0 >= 0.0) {
            if (lngp < 0.0) {
                lngp += 360.0;
            }
            else if (lngp > 360.0) {
                lngp -= 360.0;
            }
        } 
        else {
            if (lngp > 0.0) {
                lngp -= 360.0;
            }
            else if (lngp < -360.0) {
                lngp += 360.0;
            }
        }
        
    } // endelse: fiducial point was away from native pole
    
    // Store LATPOLEa
    m_ref[3] = latp;

    // Set the Euler angles
    m_euler[0] = lngp;
    m_euler[1] = 90.0 - latp;
    m_euler[2] = phip;
    sincosd(m_euler[1], &m_euler[4], &m_euler[3]);
    
    // Signal if |latitude| is preserved
    m_isolat = (m_euler[4] == 0.0);

    // Check for ill-conditioned parameters
    if (std::abs(latp) > 90.0+tol) {
        std::string message;
        message  = "abs(latp) > 90, coordinate transformation is ill-conditioned.";
        throw GException::wcs_invalid_parameter(G_CEL_SET, message);
    }

    // Celestial parameter have been set
    m_celset = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Pixel-to-world celestial transformation
 *
 * @param[in] nx X pixel vector length.
 * @param[in] ny Y pixel vector length (0=no replication, ny=nx).
 * @param[in] sxy Input vector step.
 * @param[in] sll Output vector step.
 * @param[in] x Vector of projected x coordinates.
 * @param[in] y Vector of projected y coordinates.
 * @param[out] phi Longitude of the projected point in native spherical
 *                 coordinates [deg].
 * @param[out] theta Latitude of the projected point in native spherical
 *                   coordinates [deg].
 * @param[out] lng Celestial longitude of the projected point [deg].
 * @param[out] lat Celestial latitude of the projected point [deg].
 * @param[out] stat Status return value for each vector element
 *                  (0=success, 1=invalid).
 *
 * This method transforms (x,y) coordinates in the plane of projection to
 * celestial coordinates (lng,lat). The method has been adapted from the
 * wcslib function cel.c::celx2s.
 ***************************************************************************/
void GWcslib::cel_x2s(int nx, int ny, int sxy, int sll,
                      const double* x, const double *y,
                      double* phi, double* theta,
                      double* lng, double* lat, int* stat) const
{
    // Initialize celestial transformations if required
    if (!m_celset) {
        cel_set();
    }

    // Apply spherical deprojection
    prj_x2s(nx, ny, sxy, 1, x, y, phi, theta, stat);
    
    // Compute number of Phi values
    int nphi = (ny > 0) ? (nx*ny) : nx;
    
    // Compute celestial coordinates
    sph_x2s(nphi, 0, 1, sll, phi, theta, lng, lat);

    // Return
    return;
}


/***********************************************************************//**
 * @brief World-to-pixel celestial transformation
 *
 * @param[in] nlng Longitude vector length.
 * @param[in] nlat Latitude vector length (0=no replication, nlat=nlng).
 * @param[in] sll Input vector step.
 * @param[in] sxy Output vector step.
 * @param[in] lng Celestial longitude of the projected point [deg].
 * @param[in] lat Celestial latitude of the projected point [deg].
 * @param[out] phi Longitude of the projected point in native spherical
 *                 coordinates [deg].
 * @param[out] theta Latitude of the projected point in native spherical
 *                   coordinates [deg].
 * @param[out] x Vector of projected x coordinates.
 * @param[out] y Vector of projected y coordinates.
 * @param[out] stat Status return value for each vector element
 *                  (0=success, 1=invalid).
 *
 * This method transforms (x,y) coordinates in the plane of projection to
 * celestial coordinates (lng,lat). The method has been adapted from the
 * wcslib function cel.c::celx2s.
 ***************************************************************************/
void GWcslib::cel_s2x(int nlng, int nlat, int sll, int sxy,
                      const double* lng, const double* lat,
                      double* phi, double* theta,
                      double* x, double* y, int* stat) const
{
    // Initialize celestial transformations if required
    if (!m_celset) {
        cel_set();
    }

    // Compute native coordinates
    sph_s2x(nlng, nlat, sll, 1, lng, lat, phi, theta);

    // Constant celestial latitude -> constant native latitude
    int nphi;
    int ntheta;
    if (m_isolat) {
        nphi   = nlng;
        ntheta = nlat;
    } 
    else {
        nphi   = (nlat > 0) ? (nlng*nlat) : nlng;
        ntheta = 0;
    }

    // Apply spherical projection
    prj_s2x(nphi, ntheta, 1, sxy, phi, theta, x, y, stat);
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =         Spherical transformation methods adapted from wcslib            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Rotation in the pixel-to-world direction
 *
 * @param[in] nphi Phi vector length.
 * @param[in] ntheta Theta vector length (0=no replication).
 * @param[in] spt Input vector step.
 * @param[in] sll Output vector step.
 * @param[in] phi Longitude in the native coordinate system of the
 *                projection [deg].
 * @param[in] theta Latitude in the native coordinate system of the
 *                  projection [deg].
 * @param[out] lng Celestial longitude [deg].
 * @param[out] lat Celestial latitude [deg].
 *
 * This method has been adapted from the wcslib function sph.c::sphx2s().
 * The interface follows very closely that of wcslib.
 ***************************************************************************/
void GWcslib::sph_x2s(int nphi, int ntheta, int spt, int sll,
                   const double* phi, const double* theta,
                   double* lng, double* lat) const
{
    // Set tolerance
    const double tol = 1.0e-5;
    
    // Set value replication length mphi,mtheta
    int mphi;
    int mtheta;
    if (ntheta > 0) {
        mphi   = nphi;
        mtheta = ntheta;
    }
    else {
        mphi   = 1;
        mtheta = 1;
        ntheta = nphi;
    }

    // Check for a simple change in origin of longitude
    if (m_euler[4] == 0.0) {
        
        // Initialise pointers
        double*       lngp   = lng;
        double*       latp   = lat;
        const double* phip   = phi;
        const double* thetap = theta;
            
        // Case A: ...
        if (m_euler[1] == 0.0) {
        
            // ...
            double dlng = fmod(m_euler[0] + 180.0 - m_euler[2], 360.0);

            // ...
            for (int itheta = 0; itheta < ntheta; ++itheta, phip += spt, thetap += spt) {
                for (int iphi = 0; iphi < mphi; ++iphi, lngp += sll, latp += sll) {

                    // Shift longitude
                    *lngp = *phip + dlng;
                    *latp = *thetap;

                    // Normalize the celestial longitude
                    if (m_euler[0] >= 0.0) {
                        if (*lngp < 0.0) {
                            *lngp += 360.0;
                        }
                    } 
                    else {
                        if (*lngp > 0.0) {
                            *lngp -= 360.0;
                        }
                    }
                    if (*lngp > 360.0) {
                        *lngp -= 360.0;
                    }
                    else if (*lngp < -360.0) {
                        *lngp += 360.0;
                    }
                        
                } // endfor: looped over phi
            } // endfor: looped over theta

        } // endif: ...
        
        // Case B: ...
        else {
            
            // ...
            double dlng = fmod(m_euler[0] + m_euler[2], 360.0);

            // ...
            for (int itheta = 0; itheta < ntheta; ++itheta, phip += spt, thetap += spt) {
                for (int iphi = 0; iphi < mphi; ++iphi, lngp += sll, latp += sll) {

                    // Shift longitude and inverte latitude
                    *lngp = dlng - *phip;
                    *latp = -(*thetap);

                    // Normalize the celestial longitude
                    if (m_euler[0] >= 0.0) {
                        if (*lngp < 0.0) {
                            *lngp += 360.0;
                        }
                    }
                    else {
                        if (*lngp > 0.0) {
                            *lngp -= 360.0;
                        }
                    }
                    if (*lngp > 360.0) {
                        *lngp -= 360.0;
                    }
                    else if (*lngp < -360.0) {
                        *lngp += 360.0;
                    }
                        
                } // endfor: looped over phi
            } // endfor: looped over theta
        } // endelse: ...

    } // endif: there was a simple change of longitude

    // ... otherwise to general transformation
    else {

        // Do phi dependency
        const double* phip   = phi;
        int           rowoff = 0;
        int           rowlen = nphi*sll;
        for (int iphi = 0; iphi < nphi; ++iphi, rowoff += sll, phip += spt) {
            double  dphi = *phip - m_euler[2];
            double* lngp = lng + rowoff;
            for (int itheta = 0; itheta < mtheta; ++itheta, lngp += rowlen) {
                *lngp = dphi;
            }
        }

        // Do theta dependency
        const double* thetap = theta;
        double*       lngp   = lng;
        double*       latp   = lat;
        for (int itheta = 0; itheta < ntheta; ++itheta, thetap += spt) {
            
            // ...
            double sinthe;
            double costhe;
            sincosd(*thetap, &sinthe, &costhe);
            
            // ...
            double costhe3 = costhe * m_euler[3];
            double costhe4 = costhe * m_euler[4];
            double sinthe3 = sinthe * m_euler[3];
            double sinthe4 = sinthe * m_euler[4];

            // Loop over Phi
            for (int iphi = 0; iphi < mphi; ++iphi, lngp += sll, latp += sll) {
            
                //
                double dphi = *lngp;
                double sinphi;
                double cosphi;
                sincosd(dphi, &sinphi, &cosphi);

                // Compute the celestial longitude
                double x =  sinthe4 - costhe3*cosphi;
                double y = -costhe*sinphi;
                
                // Rearrange longitude formula to reduce roundoff errors
                if (std::abs(x) < tol) {
                    x = -cosd(*thetap + m_euler[1]) + costhe3*(1.0 - cosphi);
                }

                // Compute longitude shift
                double dlng;
                if (x != 0.0 || y != 0.0) {
                    dlng = atan2d(y, x);
                }
                else {
                    dlng = (m_euler[1] < 90.0) ? dphi + 180.0 : -dphi;
                }
                
                // Set celestial longitude
                *lngp = m_euler[0] + dlng;

                // Normalize the celestial longitude
                if (m_euler[0] >= 0.0) {
                    if (*lngp < 0.0) {
                        *lngp += 360.0;
                    }
                }
                else {
                    if (*lngp > 0.0) {
                        *lngp -= 360.0;
                    }
                }
                if (*lngp > 360.0) {
                    *lngp -= 360.0;
                }
                else if (*lngp < -360.0) {
                    *lngp += 360.0;
                }

                // Compute the celestial latitude. First handle the case
                // of longitude shifts by 180 deg 
                if (fmod(dphi,180.0) == 0.0) {
                    *latp = *thetap + cosphi*m_euler[1];
                    if (*latp >  90.0) *latp =  180.0 - *latp;
                    if (*latp < -90.0) *latp = -180.0 - *latp;
                }
                
                // ... then handle the case of general longitude shifts
                else {
                    // Use alternative formulae for greater accuracy
                    double z = sinthe3 + costhe4*cosphi;
                    if (z > 0.99) {
                        *latp = acosd(sqrt(x*x+y*y));
                    }
                    else if (z < -0.99) {
                        *latp = -acosd(sqrt(x*x+y*y));
                    }
                    else {
                        *latp = asind(z);
                    }
                } // endelse: general longitude shift
                
            } // endfor: looped over phi
  
        } // endfor: looped over theta
    
    } // endelse: did general transformation

    // Return
    return;
}


/***********************************************************************//**
 * @brief Rotation in the pixel-to-world direction
 *
 * @param[in] nlng Longitude vector length.
 * @param[in] nlat Latitude vector length (0=no replication).
 * @param[in] sll Input vector step.
 * @param[in] spt Output vector step.
 * @param[in] lng Celestial longitude [deg].
 * @param[in] lat Celestial latitude [deg].
 * @param[out] phi Longitude in the native coordinate system of the
 *                 projection [deg].
 * @param[out] theta Latitude in the native coordinate system of the
 *                   projection [deg].
 *
 * This method has been adapted from the wcslib function sph.c::sphs2x().
 * The interface follows very closely that of wcslib.
 ***************************************************************************/
void GWcslib::sph_s2x(int nlng, int nlat, int sll, int spt,
                      const double* lng, const double* lat,
                      double* phi, double* theta) const
{
    // Set tolerance
    const double tol = 1.0e-5;

    // Set value replication length mlng,mlat
    int mlng;
    int mlat;
    if (nlat > 0) {
        mlng = nlng;
        mlat = nlat;
    }
    else {
        mlng = 1;
        mlat = 1;
        nlat = nlng;
    }

    // Check for a simple change in origin of longitude
    if (m_euler[4] == 0.0) {
        
        // Initialise pointers
        const double* lngp   = lng;
        const double* latp   = lat;
        double*       phip   = phi;
        double*       thetap = theta;

        // Case A: ...
        if (m_euler[1] == 0.0) {
        
            // Compute longitude shift
            double dphi = fmod(m_euler[2] - 180.0 - m_euler[0], 360.0);

            // Apply longitude shift
            for (int ilat = 0; ilat < nlat; ++ilat, lngp += sll, latp += sll) {
                for (int ilng = 0; ilng < mlng; ++ilng, phip += spt, thetap += spt) {
                
                    // Shift longitude, keep Theta
                    *phip   = fmod(*lngp + dphi, 360.0);
                    *thetap = *latp;

                    // Normalize the native longitude
                    if (*phip > 180.0) {
                        *phip -= 360.0;
                    }
                    else if (*phip < -180.0) {
                        *phip += 360.0;
                    }
        
                } // endfor: looped over longitude
            } // endfor: looped over latitude

        } // endif: Case A
        
        // Case B: ...
        else {
            
            // Compute longitude shift
            double dphi = fmod(m_euler[2] + m_euler[0], 360.0);

            // Apply longitude shift
            for (int ilat = 0; ilat < nlat; ++ilat, lngp += sll, latp += sll) {
                for (int ilng = 0; ilng < mlng; ++ilng, phip += spt, thetap += spt) {
                
                    // Shift longitude, flip Theta
                    *phip   = fmod(dphi - *lngp, 360.0);
                    *thetap = -(*latp);

                    // Normalize the native longitude
                    if (*phip > 180.0) {
                        *phip -= 360.0;
                    }
                    else if (*phip < -180.0) {
                        *phip += 360.0;
                    }
          
                } // endfor: looped over longitude
            } // endfor: looped over latitude
            
        } // endelse: Case B

    } // endif: simple change in origin
    
    // ... otherwise compute general transformation
    else {

        // Do lng dependency
        const double* lngp   = lng;
        int           rowoff = 0;
        int           rowlen = nlng * spt;
        for (int ilng = 0; ilng < nlng; ++ilng, rowoff += spt, lngp += sll) {
            double  dlng   = *lngp - m_euler[0];
            double* phip   = phi + rowoff;
            double* thetap = theta;
            for (int ilat = 0; ilat < mlat; ++ilat, phip += rowlen) {
                *phip = dlng;
            }
        }

        // Do lat dependency
        const double* latp   = lat;
        double*       phip   = phi;
        double*       thetap = theta;
        for (int ilat = 0; ilat < nlat; ++ilat, latp += sll) {
            
            // ...
            double sinlat;
            double coslat;
            sincosd(*latp, &sinlat, &coslat);
            double coslat3 = coslat*m_euler[3];
            double coslat4 = coslat*m_euler[4];
            double sinlat3 = sinlat*m_euler[3];
            double sinlat4 = sinlat*m_euler[4];

            // Loop over longitudes
            for (int ilng = 0; ilng < mlng; ++ilng, phip += spt, thetap += spt) {
            
                // ...
                double dlng = *phip;
                double sinlng;
                double coslng;
                sincosd(dlng, &sinlng, &coslng);

                // Compute the native longitude
                double x = sinlat4 - coslat3*coslng;
                double y = -coslat*sinlng;
                
                // Rearrange formula to reduce roundoff errors
                if (std::abs(x) < tol) {
                    x = -cosd(*latp+m_euler[1]) + coslat3*(1.0 - coslng);
                }

                // Compute Phi shift
                double dphi;
                if (x != 0.0 || y != 0.0) {
                    dphi = atan2d(y, x);
                } 
                else { // Change of origin of longitude
                    if (m_euler[1] < 90.0) {
                        dphi =  dlng - 180.0;
                    }
                    else {
                        dphi = -dlng;
                    }
                }
                
                // Set Phi
                *phip = fmod(m_euler[2] + dphi, 360.0);

                // Normalize the native longitude
                if (*phip > 180.0) {
                    *phip -= 360.0;
                }
                else if (*phip < -180.0) {
                    *phip += 360.0;
                }

                // Compute the native latitude. First handle the case
                // of longitude shifts by 180 deg
                if (fmod(dlng, 180.0) == 0.0) {
                    *thetap = *latp + coslng*m_euler[1];
                    if (*thetap >  90.0) *thetap =  180.0 - *thetap;
                    if (*thetap < -90.0) *thetap = -180.0 - *thetap;
                }
                
                // ... then handle the case of general longitude shifts
                else {
                    // Use alternative formulae for greater accuracy
                    double z = sinlat3 + coslat4*coslng;
                    if (z > 0.99) {
                        *thetap = acosd(sqrt(x*x+y*y));
                    }
                    else if (z < -0.99) {
                        *thetap = -acosd(sqrt(x*x+y*y));
                    }
                    else {
                        *thetap = asind(z);
                    }
                }
      
            } // endfor: looped over longitudes
    
        } // endfor: looped over latitude
  
    } // endelse: handled general case

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                   Spectral methods adapted from wcslib                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise spectral transformation parameters
 *
 * Code adapted from spc.c::spcini(). The actual version of the code does
 * nothing.
 ***************************************************************************/
void GWcslib::spc_ini(void)
{
   
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =             Linear transformation methods adapted from wcslib           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise linear transformation parameters
 *
 * Code adapted from lin.c::linini().
 ***************************************************************************/
void GWcslib::lin_ini(int naxis)
{
    // Initialise parameters
    m_linset = false;
    m_unity  = false;
    
    // Clear vectors
    m_crpix.clear();
    m_pc.clear();
    m_cdelt.clear();
    m_piximg.clear();
    m_imgpix.clear();
    
    // Continue only if there are axes
    if (naxis > 0) {
   
        // Loop over axes
        for (int i = 0; i < naxis; ++i) {
        
            // CRPIXja defaults to 0.0
            m_crpix.push_back(0.0);
            
            // PCi_ja defaults to the unit matrix
            for (int j = 0; j < naxis; ++j) {
                if (j == i) {
                    m_pc.push_back(1.0);
                }
                else {
                    m_pc.push_back(0.0);
                }
            }
            
            // CDELTia defaults to 1.0
            m_cdelt.push_back(1.0);
        
        } // endfor: looped over axes
        
    } // endif: there were axes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise linear transformation parameters
 *
 * Code adapted from lin.c::linset(). It requires m_pc to be set and to be
 * of size m_naxis*m_naxis.
 ***************************************************************************/
void GWcslib::lin_set(void) const
{
    // Clear transformation matrices
    m_piximg.clear();
    m_imgpix.clear();
    
    // Check for unity PC matrix
    m_unity = true;
    for (int i = 0, index = 0; i < m_naxis; ++i) {
        for (int j = 0; j < m_naxis; ++j, ++index) {
            if (j == i) {
                if (m_pc[index] != 1.0) {
                    m_unity = false;
                    break;
                }
            }
            else {
                if (m_pc[index] != 0.0) {
                    m_unity = false;
                    break;
                }
            }
        }
        if (!m_unity) {
            break;
        }
    }
    
    // Debug option: force PC usage for testing purposes
    #if defined(G_LIN_MATINV_FORCE_PC)
    m_unity = false;
    std::cout << "DEBUG: Force PC usage in linear transformations." << std::endl;
    #endif
    
    // Compute transformation matrices for non-uniform PC matrix
    if (!m_unity) {
    
        // Compute the pixel-to-image transformation matrix
        for (int i = 0, index = 0; i < m_naxis; ++i) {
            for (int j = 0; j < m_naxis; ++j, ++index) {
                m_piximg.push_back(m_cdelt[i] * m_pc[index]);
            }
        }
    
        // Compute the image-to-pixel transformation matrix
        lin_matinv(m_piximg, m_imgpix);
        
    }
    
    // Signal that linear transformation matrix has been set
    m_linset = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Pixel-to-world linear transformation
 *
 * @param[in] ncoord Number of coordinates.
 * @param[in] nelem Vector length of each coordinate (>=m_naxis).
 * @param[in] pixcrd Array of pixel coordinates (ncoord*nelem).
 * @param[out] imgcrd Array of intermediate world coordinates (ncoord*nelem).
 *
 * Transforms pixel coordinates to intermediate world coordinates. This
 * method has been adapted from lin.c::linp2x(). The method performs distinct
 * computations depending of whether the PC matrix is unity or not, avoiding
 * a matrix multiplication when it is not necessary.
 ***************************************************************************/
void GWcslib::lin_p2x(int ncoord, int nelem, const double* pixcrd, double* imgcrd) const
{
    // Initialize linear transformations if required
    if (!m_linset) {
        lin_set();
    }
    
    // Convert pixel coordinates to intermediate world coordinates
    const double* pix = pixcrd;
    double*       img = imgcrd;

    // Case A: we have a unity PC matrix
    if (m_unity) {
    
        // Loop over all coordinates
        for (int k = 0; k < ncoord; ++k) {
            
            // Transform the first m_naxis elements
            for (int i = 0; i < m_naxis; ++i) {
                *(img++) = m_cdelt[i] * (*(pix++) - m_crpix[i]);
            }
            
            // Go to next coordinate
            pix += (nelem - m_naxis);
            img += (nelem - m_naxis);
            
        } // endfor: looped over all coordinates

    } 
    
    // Case B: we have a non-unity PC matrix
    else {
    
        // Loop over all coordinates
        for (int k = 0; k < ncoord; ++k) {
        
            // Clear first m_naxis elements
            for (int i = 0; i < m_naxis; ++i) {
                img[i] = 0.0;
            }

            // Perform matrix multiplication (column-wise multiplication
            // allows this to be cached)
            for (int j = 0; j < m_naxis; ++j) {
                double  temp   = *(pix++) - m_crpix[j];
                for (int i = 0, ji = j; i < m_naxis; ++i, ji += m_naxis) {
                    img[i] += m_piximg[ji] * temp;
                }
            }

            // Go to next coordinate
            pix += (nelem - m_naxis);
            img += nelem;
    
        } // endfor: looped over all coordinates
        
    } // endelse: PC matrix was not unity
        
    // Return
    return;
}


/***********************************************************************//**
 * @brief World-to-pixel linear transformation
 *
 * @param[in] ncoord Number of coordinates.
 * @param[in] nelem Vector length of each coordinate (>=m_naxis).
 * @param[in] imgcrd Array of intermediate world coordinates (ncoord*nelem).
 * @param[out] pixcrd Array of pixel coordinates (ncoord*nelem).
 *
 * Transforms intermediate world coordinates to pixel coordinates. This
 * method has been adapted from lin.c::linx2p(). The method performs distinct
 * computations depending of whether the PC matrix is unity or not, avoiding
 * a matrix multiplication when it is not necessary.
 ***************************************************************************/
void GWcslib::lin_x2p(int ncoord, int nelem, const double* imgcrd, double* pixcrd) const
{
    // Initialize linear transformations if required
    if (!m_linset) {
        lin_set();
    }
    
    // Convert pixel coordinates to intermediate world coordinates
    const double* img = imgcrd;
    double*       pix = pixcrd;

    // Case A: we have a unity PC matrix
    if (m_unity) {
    
        // Loop over all coordinates
        for (int k = 0; k < ncoord; ++k) {
            
            // Transform the first m_naxis elements
            for (int i = 0; i < m_naxis; ++i) {
                *(pix++) = (*(img++) / m_cdelt[i]) + m_crpix[i];
            }
            
            // Go to next coordinate
            pix += (nelem - m_naxis);
            img += (nelem - m_naxis);
            
        } // endfor: looped over all coordinates

    } 
    
    // Case B: we have a non-unity PC matrix
    else {
    
        // Loop over all coordinates
        for (int k = 0; k < ncoord; ++k) {
        
            // Perform matrix multiplication
            for (int j = 0, ji = 0; j < m_naxis; ++j) {
                *pix = 0.0;
                for (int i = 0; i < m_naxis; ++i, ++ji) {
                    *pix += m_imgpix[ji] * img[i];
                }
                *(pix++) += m_crpix[j];
            }

            // Go to next coordinate
            pix += (nelem - m_naxis);
            img += nelem;
    
        } // endfor: looped over all coordinates
        
    } // endelse: PC matrix was not unity
        
    // Return
    return;
}


/***********************************************************************//**
 * @brief Invert linear transformation matrix
 *
 * @param[in] mat Matrix
 * @param[out] inv Inverted matrix
 *
 * @exception GException::wcs_singular_matrix
 *            Singular matrix encountered.
 *
 * Inverts a linear transformation matrix that is stored in a vector. The
 * matrix dimension is assumed to be m_naxis*m_naxis. This method has been
 * adapted from lin.c::matinv().
 ***************************************************************************/
void GWcslib::lin_matinv(const std::vector<double>& mat, std::vector<double>& inv) const
{
    // Continue only if naxis is valid
    if (m_naxis > 0) {

        // Declare internal arrays
        std::vector<int>    mxl;
        std::vector<int>    lxm;
        std::vector<double> rowmax;
        std::vector<double> lu = mat;
        
        // Clear inverse matrix
        inv.clear();
        inv.assign(m_naxis*m_naxis, 0.0);
        
        // Initialize arrays
        for (int i = 0, ij = 0; i < m_naxis; ++i) {
        
            // Vector that records row interchanges
            mxl.push_back(i);
            lxm.push_back(0);
            
            // Maximum element value in row
            rowmax.push_back(0.0);
            for (int j = 0; j < m_naxis; ++j, ++ij) {
                double dtemp = std::abs(mat[ij]);
                if (dtemp > rowmax[i]) {
                    rowmax[i] = dtemp;
                }
            }
            
            // Throw exception if matrix is singular
            if (rowmax[i] == 0.0) {
                throw GException::wcs_singular_matrix(G_LIN_MATINV, m_naxis, mat);
            }
        
        } // endfor: initialized array


        // Form the LU triangular factorization using scaled partial pivoting
        for (int k = 0; k < m_naxis; ++k) {
        
            // Decide whether to pivot
            double colmax = std::abs(lu[k*m_naxis+k]) / rowmax[k];
            int    pivot  = k;
            for (int i = k+1; i < m_naxis; ++i) {
                int    ik    = i*m_naxis + k;
                double dtemp = std::abs(lu[ik]) / rowmax[i];
                if (dtemp > colmax) {
                    colmax = dtemp;
                    pivot  = i;
                }
            }
            
            // Do we need to pivot
            if (pivot > k) {
            
                // We must pivot, interchange the rows of the design matrix
                for (int j = 0, pj = pivot*m_naxis, kj = k*m_naxis; j < m_naxis; ++j, ++pj, ++kj) {
                    double dtemp = lu[pj];
                    lu[pj]       = lu[kj];
                    lu[kj]       = dtemp;
                }

                // Amend the vector of row maxima
                double dtemp  = rowmax[pivot];
                rowmax[pivot] = rowmax[k];
                rowmax[k]     = dtemp;

                // Record the interchange for later use
                int itemp  = mxl[pivot];
                mxl[pivot] = mxl[k];
                mxl[k]     = itemp;
    
            } // endif: pivoting required

            // Gaussian elimination
            for (int i = k+1; i < m_naxis; ++i) {
                
                // Compute matrix index
                int ik = i*m_naxis + k;

                // Nothing to do if lu[ik] is zero
                if (lu[ik] != 0.0) {
                    
                    // Save the scaling factor
                    lu[ik] /= lu[k*m_naxis+k];

                    // Subtract rows
                    for (int j = k+1; j < m_naxis; ++j) {
                        lu[i*m_naxis+j] -= lu[ik]*lu[k*m_naxis+j];
                    }

                } // endif: lu[ik] was not zero
    
            } // endfor: Gaussian elimination
            
        } // endfor: formed LU triangular factorization


        // mxl[i] records which row of mat corresponds to row i of lu
        // lxm[i] records which row of lu  corresponds to row i of mat
        for (int i = 0; i < m_naxis; ++i) {
            lxm[mxl[i]] = i;
        }

        // Determine the inverse matrix
        for (int k = 0; k < m_naxis; ++k) {
        
            // ...
            inv[lxm[k]*m_naxis+k] = 1.0;

            // Forward substitution
            for (int i = lxm[k]+1; i < m_naxis; ++i) {
                for (int j = lxm[k]; j < i; ++j) {
                    inv[i*m_naxis+k] -= lu[i*m_naxis+j]*inv[j*m_naxis+k];
                }
            }

            // Backward substitution
            for (int i = m_naxis-1; i >= 0; --i) {
                for (int j = i+1; j < m_naxis; ++j) {
                    inv[i*m_naxis+k] -= lu[i*m_naxis+j]*inv[j*m_naxis+k];
                }
                inv[i*m_naxis+k] /= lu[i*m_naxis+i];
            }
            
        } // endfor: determined inverse matrix
    
    } // endif: naxis was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                  Projection methods adapted from wcslib                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise projection parameters
 *
 * Code adapted from prj.c::prjini().
 ***************************************************************************/
void GWcslib::prj_ini(void) const
{
    // Initialise parameters
    m_prjset = false;
    m_r0     = 0.0;
    m_pv[0]  = 0.0;
    m_pv[1]  = UNDEFINED;
    m_pv[2]  = UNDEFINED;
    m_pv[3]  = UNDEFINED;
    for (int k = 4; k < PVN; ++k) {
        m_pv[k] = 0.0;
    }
    m_bounds = true; // 1 in wcslib
    m_x0     = 0.0;
    m_y0     = 0.0;
    m_w.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute fiducial offset to force (x,y) = (0,0) at (phi0,theta0)
 *
 * @param[in] phi0 Fiducial longitude
 * @param[in] theta0 Fiducial latitude
 *
 * Code adapted from prj.c::prjoff().
 ***************************************************************************/
void GWcslib::prj_off(const double& phi0, const double& theta0) const
{
    // Initialise fiducial offsets
    m_x0 = 0.0;
    m_y0 = 0.0;
    
    // Set both to the projection-specific default if either undefined
    if (undefined(m_phi0) || undefined(m_theta0)) {
        m_phi0   = phi0;
        m_theta0 = theta0;
    }
    
    // ... otherwise compute (x,y) for (phi0,theta0)
    else {
        // Get (x,y) at (phi0,theta0) from projection
        int    stat;
        double x0;
        double y0;
        prj_s2x(1, 1, 1, 1, &m_phi0, &m_theta0, &x0, &y0, &stat);
        m_x0 = x0;
        m_y0 = y0;
    }
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                 Trigonometric methods adapted from wcslib               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Compute cosine of angle in degrees
 *
 * @param[in] angle Angle in degrees
 *
 * Code adapted from wcstrig.c::cosd().
 ***************************************************************************/
double GWcslib::cosd(const double& angle) const
{
    // Check for rounding errors
    if (fmod(angle, 90.0) == 0.0) {
        int i = std::abs((int)std::floor(angle/90.0 + 0.5))%4;
        switch (i) {
        case 0:
            return 1.0;
        case 1:
            return 0.0;
        case 2:
            return -1.0;
        case 3:
            return 0.0;
        }
    }

    // Return cosine
    return std::cos(angle*deg2rad);
}


/***********************************************************************//**
 * @brief Compute sine of angle in degrees
 *
 * @param[in] angle Angle in degrees
 *
 * Code adapted from wcstrig.c::sind().
 ***************************************************************************/
double GWcslib::sind(const double& angle) const
{
    // Check for rounding errors
    if (fmod(angle, 90.0) == 0.0) {
        int i = std::abs((int)std::floor(angle/90.0 - 0.5))%4;
        switch (i) {
        case 0:
            return 1.0;
        case 1:
            return 0.0;
        case 2:
            return -1.0;
        case 3:
            return 0.0;
        }
    }

    // Return sine
    return std::sin(angle*deg2rad);
}


/***********************************************************************//**
 * @brief Compute tangens of angle in degrees
 *
 * @param[in] angle Angle in degrees
 *
 * Code adapted from wcstrig.c::tand().
 ***************************************************************************/
double GWcslib::tand(const double& angle) const
{
    // Check for rounding errors
    double resid = fmod(angle, 360.0);
    if (resid == 0.0 || std::abs(resid) == 180.0) {
        return 0.0;
    }
    else if (resid == 45.0 || resid == 225.0) {
        return 1.0;
    }
    else if (resid == -135.0 || resid == -315.0) {
        return -1.0;
    }

    // Return tangens
    return std::tan(angle*deg2rad);
}


/***********************************************************************//**
 * @brief Compute arc cosine in degrees
 *
 * @param[in] value Value
 *
 * Code adapted from wcstrig.c::acosd().
 ***************************************************************************/
double GWcslib::acosd(const double& value) const
{
    // Domain tolerance
    const double wcstrig_tol = 1.0e-10;
    
    // Check for rounding errors
    if (value >= 1.0) {
        if (value-1.0 <  wcstrig_tol) {
            return 0.0;
        }
    } 
    else if (value == 0.0) {
        return 90.0;
    }
    else if (value <= -1.0) {
        if (value+1.0 > -wcstrig_tol) {
            return 180.0;
        }
    }

    // Return arc cosine
    return std::acos(value)*rad2deg;
}


/***********************************************************************//**
 * @brief Compute arc sine in degrees
 *
 * @param[in] value Value
 *
 * Code adapted from wcstrig.c::asind().
 ***************************************************************************/
double GWcslib::asind(const double& value) const
{
    // Domain tolerance
    const double wcstrig_tol = 1.0e-10;
    
    // Check for rounding errors
    if (value <= -1.0) {
        if (value+1.0 > -wcstrig_tol) {
            return -90.0;
        }
    } 
    else if (value == 0.0) {
        return 0.0;
    }
    else if (value >= 1.0) {
        if (value-1.0 <  wcstrig_tol) {
            return 90.0;
        }
    }

    // Return arc sine
    return std::asin(value)*rad2deg;
}


/***********************************************************************//**
 * @brief Compute arc tangens in degrees
 *
 * @param[in] value Value
 *
 * Code adapted from wcstrig.c::atand().
 ***************************************************************************/
double GWcslib::atand(const double& value) const
{
    // Check for rounding errors
    if (value == -1.0) {
        return -45.0;
    }
    else if (value == 0.0) {
        return 0.0;
    }
    else if (value == 1.0) {
        return 45.0;
    }

    // Return arc sine
    return std::atan(value)*rad2deg;
}


/***********************************************************************//**
 * @brief Compute arc tangens in degrees
 *
 * @param[in] y Nominator
 * @param[in] x Denominator
 *
 * Code adapted from wcstrig.c::atan2d().
 ***************************************************************************/
double GWcslib::atan2d(const double& y, const double& x) const
{
    // Check for rounding errors
    if (y == 0.0) {
        if (x >= 0.0) {
            return 0.0;
        }
        else if (x < 0.0) {
            return 180.0;
        }
    }
    else if (x == 0.0) {
        if (y > 0.0) {
            return 90.0;
        }
        else if (y < 0.0) {
            return -90.0;
        }
    }

    // Return arc sine
    return std::atan2(y,x)*rad2deg;
}


/***********************************************************************//**
 * @brief Compute sine and cosine of angle in degrees
 *
 * @param[in] angle Angle [degrees].
 * @param[out] s Sine of angle.
 * @param[out] c Cosine of angle.
 *
 * Code adapted from wcstrig.c::sincosd().
 ***************************************************************************/
void GWcslib::sincosd(const double& angle, double *s, double *c) const

{
    // Check for rounding errors
    if (fmod(angle, 90.0) == 0.0) {
        int i = std::abs((int)std::floor(angle/90.0 + 0.5))%4;
        switch (i) {
        case 0:
            *s = 0.0;
            *c = 1.0;
            return;
        case 1:
            *s = (angle > 0.0) ? 1.0 : -1.0;
            *c = 0.0;
            return;
        case 2:
            *s =  0.0;
            *c = -1.0;
            return;
        case 3:
            *s = (angle > 0.0) ? -1.0 : 1.0;
            *c = 0.0;
            return;
        }
    }
  
    // Compute sine and cosine
    *s = std::sin(angle*deg2rad);
    *c = std::cos(angle*deg2rad);
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
