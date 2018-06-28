/***************************************************************************
 *               GWcsAIT.cpp - Aitoff (AIT) projection class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Juergen Knoedlseder                         *
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
 * @file GWcsAIT.cpp
 * @brief Aitoff (AIT) projection class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GTools.hpp"
#include "GWcsAIT.hpp"
#include "GWcsRegistry.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PRJ_X2S    "GWcsAIT::prj_x2s(int, int, int, int, double*, double*,"\
                                                   " double*, double*, int*)"
#define G_PRJ_S2X    "GWcsAIT::prj_s2x(int, int, int, int, double*, double*,"\
                                                   " double*, double*, int*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_PRJ                             //!< Debug GWcsAIT::prj_x2s

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GWcsAIT      g_wcs_ait_seed;
const GWcsRegistry g_wcs_ait_registry(&g_wcs_ait_seed);


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcsAIT::GWcsAIT(void) : GWcs()
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
GWcsAIT::GWcsAIT(const std::string& coords,
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
GWcsAIT::GWcsAIT(const GWcsAIT& wcs) : GWcs(wcs)
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
GWcsAIT::~GWcsAIT(void)
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
GWcsAIT& GWcsAIT::operator=(const GWcsAIT& wcs)
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
void GWcsAIT::clear(void)
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
 * @return Pointer to deep copy of Aitoff projection.
 ***************************************************************************/
GWcsAIT* GWcsAIT::clone(void) const
{
    return new GWcsAIT(*this);
}



/***********************************************************************//**
 * @brief Print WCS information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing WCS information.
 ***************************************************************************/
std::string GWcsAIT::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GWcsAIT ===");

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
void GWcsAIT::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
void GWcsAIT::copy_members(const GWcsAIT& wcs)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsAIT::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup of projection
 *
 * This method sets up the projection information. The method has been
 * adapted from the wcslib function prj.c::aitset.
 *
 * @exception GException::wcs_invalid_parameter
 *            PV(1) or PV(2) are invalid.
 *
 *   Given and/or returned:
 *      m_r0      Reset to 180/pi if 0.
 *      m_phi0    Reset to  0.0 if undefined.
 *      m_theta0  Reset to  0.0 if undefined.
 *
 *   Returned:
 *      m_x0      Fiducial offset in x.
 *      m_y0      Fiducial offset in y.
 *      m_w[0]    2*r0**2
 *      m_w[1]    1/(2*r0)**2
 *      m_w[2]    1/(4*r0)**2
 *      m_w[3]    1/(2*r0)
 ***************************************************************************/
void GWcsAIT::prj_set(void) const
{
    // Signal that projection has been set (needs to be done before calling
    // the prj_off() method to avoid an endless loop)
    m_prjset = true;

    // Initialise projection parameters
    m_w.assign(4, 0.0);
    
    // Set undefined parameters
    if (m_r0 == 0.0) {
        m_r0 = gammalib::rad2deg;
    }
    
    // Precompute 
    m_w[0] = 2.0 * m_r0 * m_r0;
    m_w[1] = 1.0 / (2.0 * m_w[0]);
    m_w[2] = m_w[1] / 4.0;
    m_w[3] = 1.0 / (2.0 * m_r0);
    
    // Compute fiducial offset
    prj_off(0.0, 0.0);
    
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
 * This method has been adapted from the wcslib function prj.c::aitx2s().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previously (as this will be done by the constructor).
 ***************************************************************************/
void GWcsAIT::prj_x2s(int nx, int ny, int sxy, int spt, 
                      const double* x, const double* y,
                      double* phi, double* theta, int* stat) const
{
    // Set tolerance
    const double tol = 1.0e-13;

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

    // Initialise status code and statistics
    int status    = 0;
    int n_invalid = 0;
    
    // Do x dependence
    const double* xp     = x;
    int           rowoff = 0;
    int           rowlen = nx * spt;
    for (int ix = 0; ix < nx; ++ix, rowoff += spt, xp += sxy) {
        double  xj     = *xp + m_x0;
        double  s      = 1.0 - xj * xj * m_w[2];
        double  t      = xj * m_w[3];
        double* phip   = phi   + rowoff;
        double* thetap = theta + rowoff;
        for (int iy = 0; iy < my; ++iy, phip += rowlen, thetap += rowlen) {
            *phip   = s;
            *thetap = t;
        }
    }

    // Do y dependence
    const double* yp     = y;
    double*       phip   = phi;
    double*       thetap = theta;
    int*          statp  = stat;
    for (int iy = 0; iy < ny; ++iy, yp += sxy) {
        double yj  = *yp + m_y0;
        double yj2 = yj * yj * m_w[1];
        for (int ix = 0; ix < mx; ++ix, phip += spt, thetap += spt) {

            // Initialise status
            int istat = 0;

            // Compute Phi
            double s = *phip - yj2;
            if (s < 0.5) {
                if (s < 0.5-tol) {
                    istat  = 1;
                    status = 3;
                    n_invalid++;
                    #if defined(G_DEBUG_PRJ)
                    std::cout << "prj_x2s(Phi)..:";
                    std::cout << " nx=" << nx;
                    std::cout << " ny=" << ny;
                    std::cout << " ix=" << ix;
                    std::cout << " iy=" << iy;
                    std::cout << " yp=" << *yp;
                    std::cout << " yj=" << yj;
                    std::cout << " yj2=" << yj2;
                    std::cout << " phip=" << *phip;
                    std::cout << " s=" << s;
                    std::cout << std::endl;
                    #endif
                }
                s = 0.5;
            }
            double z  = std::sqrt(s);
            double x0 = 2.0*z*z - 1.0;
            double y0 = z*(*thetap);
            if (x0 == 0.0 && y0 == 0.0) {
                *phip = 0.0;
            }
            else {
                *phip = 2.0 * gammalib::atan2d(y0, x0);
            }

            // Compute Theta
            double t = z*yj/m_r0;
            if (std::abs(t) > 1.0) {
                if (std::abs(t) > 1.0+tol) {
                    if (istat == 0) {
                        istat  = 1;
                        status = 3;
                        n_invalid++;
                    }
                    #if defined(G_DEBUG_PRJ)
                    std::cout << "prj_x2s(Theta):";
                    std::cout << " x0=" << x0;
                    std::cout << " y0=" << y0;
                    std::cout << " phip=" << *phip;
                    std::cout << " z=" << z;
                    std::cout << " yj=" << yj;
                    std::cout << " m_r0=" << m_r0;
                    std::cout << " t=" << t;
                    std::cout << std::endl;
                    #endif
                }
                t = (t < 0) ? -90.0 : 90.0;
            }
            else {
                t = gammalib::asind(t);
            }
            *thetap    = t;
            *(statp++) = istat;
        }
    }

    // Handle status code
    if (status == 3) {
        throw GException::wcs_invalid_x_y(G_PRJ_X2S, n_invalid);
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
 * This method has been adapted from the wcslib function prj.c::aits2x().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previously (as this will be done by the constructor).
 ***************************************************************************/
void GWcsAIT::prj_s2x(int nphi, int ntheta, int spt, int sxy,
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
    int status    = 0;
    int n_invalid = 0;
    
    // Do phi dependence
    const double* phip   = phi;
    int           rowoff = 0;
    int           rowlen = nphi * sxy;
    for (int iphi = 0; iphi < nphi; ++iphi, rowoff += sxy, phip += spt) {
        double w = (*phip)/2.0;
        double sinphi;
        double cosphi;
        gammalib::sincosd(w, &sinphi, &cosphi);
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
    
        // Compute sin(theta) and cos(theta)
        double sinthe;
        double costhe;
        gammalib::sincosd(*thetap, &sinthe, &costhe);

        // Do phi dependence
        for (int iphi = 0; iphi < mphi; ++iphi, xp += sxy, yp += sxy) {
            double w   = std::sqrt(m_w[0]/(1.0 + costhe * (*yp)));
            *xp        = 2.0 * w * costhe * (*xp) - m_x0;
            *yp        = w * sinthe - m_y0;
            *(statp++) = 0;
        } // endfor: phi

    } // endfor: theta
  
    // Handle status code
    if (status == 4) {
        throw GException::wcs_invalid_phi_theta(G_PRJ_S2X, n_invalid);
    }
    
    // Return
    return;
}
