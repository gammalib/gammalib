/***************************************************************************
 *   GWcsARC.cpp - Zenithal/azimuthal equidistant (ARC) projection class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GWcsARC.cpp
 * @brief Zenithal/azimuthal equidistant (ARC) projection class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GTools.hpp"
#include "GWcsARC.hpp"
#include "GWcsRegistry.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PRJ_X2S    "GWcsARC::prj_x2s(int, int, int, int, double*, double*,"\
                                                   " double*, double*, int*)"
#define G_PRJ_S2X    "GWcsARC::prj_s2x(int, int, int, int, double*, double*,"\
                                                   " double*, double*, int*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GWcsARC      g_wcs_arc_seed;
const GWcsRegistry g_wcs_arc_registry(&g_wcs_arc_seed);


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcsARC::GWcsARC(void) : GWcs()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel (starting from 1).
 * @param[in] crpix2 Y index of reference pixel (starting from 1).
 * @param[in] cdelt1 Increment in x direction at reference pixel [deg].
 * @param[in] cdelt2 Increment in y direction at reference pixel [deg].
 ***************************************************************************/
GWcsARC::GWcsARC(const std::string& coords,
                 const double& crval1, const double& crval2,
                 const double& crpix1, const double& crpix2,
                 const double& cdelt1, const double& cdelt2) :
                 GWcs(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2)

{
    // Initialise class members
    init_members();

    // Setup WCS derived parameters
    wcs_set();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
GWcsARC::GWcsARC(const GWcsARC& wcs) : GWcs(wcs)
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
GWcsARC::~GWcsARC(void)
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
 * @return World Coordinate System.
 ***************************************************************************/
GWcsARC& GWcsARC::operator=(const GWcsARC& wcs)
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
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GWcsARC::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GWcs::free_members();
    this->GSkyProjection::free_members();

    // Initialise members
    this->GSkyProjection::init_members();
    this->GWcs::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of Zenithal/azimuthal equidistant projection.
 ***************************************************************************/
GWcsARC* GWcsARC::clone(void) const
{
    return new GWcsARC(*this);
}



/***********************************************************************//**
 * @brief Print WCS information
 *
 * @param[in] chatter Chattiness.
 * @return String containing WCS information.
 ***************************************************************************/
std::string GWcsARC::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GWcsARC ===");

        // Append information
        result.append(wcs_print(chatter));

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
void GWcsARC::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
void GWcsARC::copy_members(const GWcsARC& wcs)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsARC::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup of projection
 *
 * This method sets up the projection information. The method has been
 * adapted from the wcslib function prj.c::arcset (version 7.2).
 *
 *   Given and/or returned:
 *      m_r0      Reset to 180/pi if 0.
 *      m_phi0    Reset to  0.0 if undefined.
 *      m_theta0  Reset to 90.0 if undefined.
 *
 *   Returned:
 *      m_x0      Fiducial offset in x.
 *      m_y0      Fiducial offset in y.
 *      m_w[0]    r0*(pi/180)
 *      m_w[1]    (180/pi)/r0
 ***************************************************************************/
void GWcsARC::prj_set(void) const
{
    // Signal that projection has been set (needs to be done before calling
    // the prj_off() method to avoid an endless loop)
    m_prjset = true;

    // Initialise projection parameters
    m_w.assign(2, 0.0);
    
    // Set undefined parameters
    if (m_r0 == 0.0) {
        m_r0   = gammalib::rad2deg;
        m_w[0] = 1.0;
        m_w[1] = 1.0;
    }
    else {
        m_w[0] = m_r0 * gammalib::deg2rad;
        m_w[1] = 1.0/m_w[0];
    }
    
    // Compute fiducial offset
    prj_off(0.0, 90.0);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Pixel-to-spherical deprojection
 *
 * @param[in] nx X vector length.
 * @param[in] ny Y vector length (0=no replication).
 * @param[in] sxy Input vector step.
 * @param[in] spt Output vector step.
 * @param[in] x Vector of projected x coordinates.
 * @param[in] y Vector of projected y coordinates.
 * @param[out] phi Longitude of the projected point in native spherical
 *                 coordinates [deg].
 * @param[out] theta Latitude of the projected point in native spherical
 *                   coordinates [deg].
 * @param[out] stat Status return value for each vector element (always 0)
 *
 * @exception GException::wcs_invalid_x_y
 *            One or more of the (x,y) coordinates were invalid, as indicated
 *            by the stat vector.
 *
 * Deproject pixel (x,y) coordinates in the plane of projection to native
 * spherical coordinates (phi,theta).
 *
 * This method has been adapted from the wcslib function prj.c::arcx2s().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previously (as this will be done by the constructor).
 ***************************************************************************/
void GWcsARC::prj_x2s(int nx, int ny, int sxy, int spt, 
                      const double* x, const double* y,
                      double* phi, double* theta, int* stat) const
{
    // Initialize projection if required
    if (!m_prjset) {
        prj_set();
    }

    // Set value replication length mx,my
    int mx;
    int my;
    if (ny > 0) {
        mx = nx;
        my = ny;
    } 
    else {
        mx = 1;
        my = 1;
        ny = nx;
    }

    // Do x dependence
    const double* xp     = x;
    int           rowoff = 0;
    int           rowlen = nx * spt;
    for (int ix = 0; ix < nx; ++ix, rowoff += spt, xp += sxy) {
        double  xj     = *xp + m_x0;
        double* phip   = phi + rowoff;
        for (int iy = 0; iy < my; ++iy, phip += rowlen) {
            *phip   = xj;
        }
    }

    // Do y dependence
    const double* yp     = y;
    double*       phip   = phi;
    double*       thetap = theta;
    int*          statp  = stat;
    for (int iy = 0; iy < ny; ++iy, yp += sxy) {
        double yj  = *yp + m_y0;
        double yj2 = yj * yj;
        for (int ix = 0; ix < mx; ++ix, phip += spt, thetap += spt) {
            double xj = *phip;
            double r  = std::sqrt(xj*xj + yj2);
            if (r == 0.0) {
                *phip   =  0.0;
                *thetap = 90.0;
            }
            else {
                *phip   = gammalib::atan2d(xj, -yj);
                *thetap = 90.0 - r * m_w[1];
            }
            *(statp++) = 0;
        }
    }

    // Do boundary checking
    int status = prj_bchk(1.0e-13, nx, my, spt, phi, theta, stat);
    
    // Handle status code
    if (status != 0) {
        std::string msg = "One or more of the (x, y) coordinates were invalid "
                          "for the "+name()+" ("+code()+") projection.";
        throw GException::invalid_value(G_PRJ_X2S, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generic spherical-to-pixel projection
 *
 * @param[in] nphi Longitude vector length.
 * @param[in] ntheta Latitude vector length (0=no replication).
 * @param[in] spt Input vector step.
 * @param[in] sxy Output vector step.
 * @param[in] phi Longitude vector of the projected point in native spherical
 *                coordinates [deg].
 * @param[in] theta Latitude vector of the projected point in native spherical
 *                  coordinates [deg].
 * @param[out] x Vector of projected x coordinates.
 * @param[out] y Vector of projected y coordinates.
 * @param[out] stat Status return value for each vector element (always 0)
 *
 * @exception GException::wcs_invalid_phi_theta
 *            One or more of the (phi,theta) coordinates were invalid, as
 *            indicated by the stat vector.
 *
 * Project native spherical coordinates (phi,theta) to pixel (x,y)
 * coordinates in the plane of projection.
 *
 * This method has been adapted from the wcslib function prj.c::arcs2x().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previously (as this will be done by the constructor).
 ***************************************************************************/
void GWcsARC::prj_s2x(int nphi, int ntheta, int spt, int sxy,
                      const double* phi, const double* theta,
                      double* x, double* y, int* stat) const
{
    // Initialize projection if required
    if (!m_prjset) {
        prj_set();
    }

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

    // Initialise status code and statistics
    //int status    = 0;
    //int n_invalid = 0;
    
    // Do phi dependence
    const double* phip   = phi;
    int           rowoff = 0;
    int           rowlen = nphi * sxy;
    for (int iphi = 0; iphi < nphi; ++iphi, rowoff += sxy, phip += spt) {
        double sinphi;
        double cosphi;
        gammalib::sincosd(*phip, &sinphi, &cosphi);
        double* xp = x + rowoff;
        double* yp = y + rowoff;
        for (int itheta = 0; itheta < mtheta; ++itheta) {
            *xp = sinphi;
            *yp = cosphi;
            xp += rowlen;
            yp += rowlen;
        }
    }

    // Do theta dependence
    const double* thetap = theta;
    double*       xp     = x;
    double*       yp     = y;
    int*          statp  = stat;
    for (int itheta = 0; itheta < ntheta; ++itheta, thetap += spt) {
    
        // Compute r
        double r = m_w[0] * (90.0 - *thetap);

        // Do phi dependence
        for (int iphi = 0; iphi < mphi; ++iphi, xp += sxy, yp += sxy) {
            *xp        =  r * (*xp) - m_x0;
            *yp        = -r * (*yp) - m_y0;
            *(statp++) = 0;
        } // endfor: phi

    } // endfor: theta

    // Handle status code
    /*
    if (status == 4) {
        throw GException::wcs_invalid_phi_theta(G_PRJ_S2X, n_invalid);
    }
    */
    
    // Return
    return;
}
