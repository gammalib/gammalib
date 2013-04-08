/***************************************************************************
 *   GWcsAZP.cpp - Zenithal/azimuthal perspective (AZP) projection class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GWcsAZP.cpp
 * @brief Zenithal/azimuthal perspective (AZP) projection class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcsAZP.hpp"
#include "GWcsRegistry.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PRJ_SET                                        "GWcsAZP::prj_set()"
#define G_PRJ_X2S "GWcsAZP::prj_x2s(int,int,int,int,double*,double*,double*,"\
                                                              "double*,int*)"
#define G_PRJ_S2X "GWcsAZP::prj_s2x(int,int,int,int,double*,double*,double*,"\
                                                              "double*,int*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GWcsAZP      g_wcs_azp_seed;
const GWcsRegistry g_wcs_azp_registry(&g_wcs_azp_seed);


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcsAZP::GWcsAZP(void) : GWcslib()
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
GWcsAZP::GWcsAZP(const std::string& coords,
                 const double& crval1, const double& crval2,
                 const double& crpix1, const double& crpix2,
                 const double& cdelt1, const double& cdelt2) :
                 GWcslib(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2)

{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
GWcsAZP::GWcsAZP(const GWcsAZP& wcs) : GWcslib(wcs)
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
GWcsAZP::~GWcsAZP(void)
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
GWcsAZP& GWcsAZP::operator= (const GWcsAZP& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

        // Copy base class members
        this->GWcslib::operator=(wcs);

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
void GWcsAZP::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GWcslib::free_members();
    this->GWcs::free_members();

    // Initialise members
    this->GWcs::init_members();
    this->GWcslib::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GWcsAZP* GWcsAZP::clone(void) const
{
    return new GWcsAZP(*this);
}



/***********************************************************************//**
 * @brief Print WCS information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing WCS information.
 ***************************************************************************/
std::string GWcsAZP::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GWcsAZP ===");

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
 *
 * This method sets up the World Coordinate System by calling wcs_set().
 ***************************************************************************/
void GWcsAZP::init_members(void)
{
    // Setup World Coordinate System
    wcs_set();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
void GWcsAZP::copy_members(const GWcsAZP& wcs)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsAZP::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup of projection
 *
 * This method sets up the projection information. The method has been
 * adapted from the wcslib function prj.c::carset.
 *
 * @exception GException::wcs_invalid_parameter
 *            PV(1) or PV(2) are invalid.
 *
 *   Given:
 *      m_pv[1]   Distance parameter, mu in units of r0.
 *      m_pv[2]   Tilt angle, gamma in degrees.
 *
 *   Given and/or returned:
 *      m_r0      Reset to 180/pi if 0.
 *      m_phi0    Reset to  0.0 if undefined.
 *      m_theta0  Reset to 90.0 if undefined.
 *
 *   Returned:
 *      m_x0      Fiducial offset in x.
 *      m_y0      Fiducial offset in y.
 *      m_w[0]    r0*(mu+1)
 *      m_w[1]    tan(gamma)
 *      m_w[2]    sec(gamma)
 *      m_w[3]    cos(gamma)
 *      m_w[4]    sin(gamma)
 *      m_w[5]    asin(-1/mu) for |mu| >= 1, -90 otherwise
 *      m_w[6]    mu*cos(gamma)
 *      m_w[7]    1 if |mu*cos(gamma)| < 1, 0 otherwise
 ***************************************************************************/
void GWcsAZP::prj_set(void) const
{
    // Initialise projection parameters
    m_w.assign(8,0.0);
    
    // Set undefined parameters
    if (undefined(m_pv[1])) m_pv[1] = 0.0;
    if (undefined(m_pv[2])) m_pv[2] = 0.0;
    if (m_r0 == 0.0)        m_r0    = rad2deg;
    
    // Precompute 
    m_w[0] = m_r0*(m_pv[1] + 1.0);
    if (m_w[0] == 0.0) {
        std::string message = "PV(1)=-1 encountered.";
        throw GException::wcs_invalid_parameter(G_PRJ_SET, message);
    }
    m_w[3] = cosd(m_pv[2]);
    if (m_w[3] == 0.0) {
        std::string message = "cos(PV(2))=0 encountered.";
        throw GException::wcs_invalid_parameter(G_PRJ_SET, message);
    }
    m_w[2] = 1.0/m_w[3];
    m_w[4] = sind(m_pv[2]);
    m_w[1] = m_w[4] / m_w[3];
    if (std::abs(m_pv[1]) > 1.0)
        m_w[5] = asind(-1.0/m_pv[1]);
    else
        m_w[5] = -90.0;
    m_w[6] = m_pv[1] * m_w[3];
    m_w[7] = (std::abs(m_w[6]) < 1.0) ? 1.0 : 0.0;
    
    // Compute fiducial offset
    prj_off(0.0, 90.0);
    
    // Signal that projection has been set
    m_prjset = true;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Cartesian-to-spherical deprojection
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
 * Deproject Cartesian (x,y) coordinates in the plane of projection to native
 * spherical coordinates (phi,theta).
 *
 * This method has been adapted from the wcslib function prj.c::carx2s().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previsouly (as this will be done by the constructor).
 ***************************************************************************/
void GWcsAZP::prj_x2s(int nx, int ny, int sxy, int spt, 
                      const double* x, const double* y,
                      double* phi, double* theta, int* stat) const
{
    // Set tolerance
    const double tol = 1.0e-13;

    // Initialize projection if required
    if (!m_prjset)
        prj_set();

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
        double  xj   = *xp + m_x0;
        double* phip = phi + rowoff;
        for (int iy = 0; iy < my; ++iy, phip += rowlen)
            *phip = xj;
    }

    // Do y dependence
    const double* yp     = y;
    double*       phip   = phi;
    double*       thetap = theta;
    int*          statp  = stat;
    for (int iy = 0; iy < ny; ++iy, yp += sxy) {
        double yj  = *yp + m_y0;
        double yc  = yj * m_w[3];
        double yc2 = yc*yc;
        double q   = m_w[0] + yj*m_w[4];
        for (int ix = 0; ix < mx; ++ix, phip += spt, thetap += spt) {
            double xj = *phip;
            double r  = sqrt(xj*xj + yc2);
            if (r == 0.0) {
                *phip      =  0.0;
                *thetap    = 90.0;
                *(statp++) =    0;

            } 
            else {
                *phip    = atan2d(xj, -yc);
                double s = r / q;
                double t = s * m_pv[1]/sqrt(s*s + 1.0);
                s = atan2d(1.0, s);
                if (std::abs(t) > 1.0) {
                    if (std::abs(t) > 1.0+tol) {
                        *thetap    = 0.0;
                        *(statp++) =   1;
                        status     =   3;
                        n_invalid++;
                        continue;
                    }
                    t = (t < 0) ? -90.0 : 90.0;
                } 
                else
                    t = asind(t);
                double a = s - t;
                double b = s + t + 180.0;
                if (a > 90.0) a -= 360.0;
                if (b > 90.0) b -= 360.0;
                *thetap    = (a > b) ? a : b;
                *(statp++) = 0;
            }
        }
    }

    // Handle status code
    if (status == 3)
        throw GException::wcs_invalid_x_y(G_PRJ_X2S, n_invalid);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Generic spherical-to-Cartesian projection
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
 * Project native spherical coordinates (phi,theta) to Cartesian (x,y)
 * coordinates in the plane of projection.
 *
 * This method has been adapted from the wcslib function prj.c::cars2x().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previsouly (as this will be done by the constructor).
 ***************************************************************************/
void GWcsAZP::prj_s2x(int nphi, int ntheta, int spt, int sxy,
                      const double* phi, const double* theta,
                      double* x, double* y, int* stat) const
{
    // Initialize projection if required
    if (!m_prjset)
        prj_set();

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
        double sinphi;
        double cosphi;
        sincosd(*phip, &sinphi, &cosphi);
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
        sincosd(*thetap, &sinthe, &costhe);

        // ...
        for (int iphi = 0; iphi < mphi; ++iphi, xp += sxy, yp += sxy) {
        
            double s = m_w[1]*(*yp);
            double t = (m_pv[1] + sinthe) + costhe*s;

            // ...
            if (t == 0.0) {
                *xp        = 0.0;
                *yp        = 0.0;
                *(statp++) = 1;
                status     = 4;
                n_invalid++;
            } 
            
            // ...
            else {
                double r = m_w[0]*costhe/t;

                // Bounds checking
                int istat = 0;
                if (m_bounds) {
                
                    // Is there overlap?
                    if (*thetap < m_w[5]) {
                        istat  = 1;
                        status = 4;
                        n_invalid++;
                    }
                    
                    // Is there divergence?
                    else if (m_w[7] > 0.0) {
                        t = m_pv[1] / sqrt(1.0 + s*s);
                        if (std::abs(t) <= 1.0) {
                            s = atand(-s);
                            t = asind(t);
                            double a = s - t;
                            double b = s + t + 180.0;
                            if (a > 90.0) a -= 360.0;
                            if (b > 90.0) b -= 360.0;
                            if (*thetap < ((a > b) ? a : b)) {
                                istat  = 1;
                                status = 4;
                                n_invalid++;
                            }
                        }
                    }
                }
                
                // Set value
                *xp        =  r*(*xp) - m_x0;
                *yp        = -r*(*yp) * m_w[2] - m_y0;
                *(statp++) = istat;
            } // endelse
        } // endfor: phi
    } // endfor: theta
  
    // Handle status code
    if (status == 4)
        throw GException::wcs_invalid_phi_theta(G_PRJ_S2X, n_invalid);
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
