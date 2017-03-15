/***************************************************************************
 *        GWcsSIN.cpp - Orthographic/synthesis (SIN) projection class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file GWcsSIN.cpp
 * @brief Orthographic/synthesis (SIN) projection class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GWcsSIN.hpp"
#include "GWcsRegistry.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PRJ_X2S    "GWcsSIN::prj_x2s(int, int, int, int, double*, double*,"\
                                                   " double*, double*, int*)"
#define G_PRJ_S2X    "GWcsSIN::prj_s2x(int, int, int, int, double*, double*,"\
                                                   " double*, double*, int*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GWcsSIN      g_wcs_sin_seed;
const GWcsRegistry g_wcs_sin_registry(&g_wcs_sin_seed);


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GWcsSIN::GWcsSIN(void) : GWcs()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Projection constructor
 *
 * @param[in] coords Coordinate system.
 * @param[in] crval1 X value of reference pixel.
 * @param[in] crval2 Y value of reference pixel.
 * @param[in] crpix1 X index of reference pixel (starting from 1).
 * @param[in] crpix2 Y index of reference pixel (starting from 1).
 * @param[in] cdelt1 Increment in x direction at reference pixel [deg].
 * @param[in] cdelt2 Increment in y direction at reference pixel [deg].
 ***************************************************************************/
GWcsSIN::GWcsSIN(const std::string& coords,
                 const double& crval1, const double& crval2,
                 const double& crpix1, const double& crpix2,
                 const double& cdelt1, const double& cdelt2) :
                 GWcs(coords, crval1, crval2, crpix1, crpix2, cdelt1, cdelt2)

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
GWcsSIN::GWcsSIN(const GWcsSIN& wcs) : GWcs(wcs)
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
GWcsSIN::~GWcsSIN(void)
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
GWcsSIN& GWcsSIN::operator=(const GWcsSIN& wcs)
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
 * @brief Clear projection
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GWcsSIN::clear(void)
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
 * @brief Clone projection
 *
 * @return Pointer to deep copy of projection.
 ***************************************************************************/
GWcsSIN* GWcsSIN::clone(void) const
{
    return new GWcsSIN(*this);
}


/***********************************************************************//**
 * @brief Print projection information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing projection information.
 ***************************************************************************/
std::string GWcsSIN::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GWcsSIN ===");

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
void GWcsSIN::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
void GWcsSIN::copy_members(const GWcsSIN& wcs)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsSIN::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup of projection
 *
 * This method sets up the projection information. The method has been
 * adapted from the wcslib function prj.c::sinset.
 *
 *   Given:
 *      m_pv[1]   Obliqueness parameter xi.
 *      m_pv[2]   Obliqueness parameter eta.
 *
 *   Given and/or returned:
 *      m_r0      Reset to 180/pi if 0.
 *      m_phi0    Reset to  0.0 if undefined.
 *      m_theta0  Reset to 90.0 if undefined.
 *
 *   Returned:
 *      m_x0      Fiducial offset in x.
 *      m_y0      Fiducial offset in y.
 *      m_w[0]    1/r0
 *      m_w[1]    xi**2 + eta**2
 *      m_w[2]    xi**2 + eta**2 + 1
 *      m_w[3]    xi**2 + eta**2 - 1
 ***************************************************************************/
void GWcsSIN::prj_set(void) const
{
    // Signal that projection has been set (needs to be done before calling
    // the prj_off() method to avoid an endless loop)
    m_prjset = true;

    // Initialise projection parameters
    m_w.assign(4,0.0);
    
    // Set undefined parameters
    if (undefined(m_pv[1])) {
        m_pv[1] = 0.0;
    }
    if (undefined(m_pv[2])) {
        m_pv[2] = 0.0;
    }
    if (m_r0 == 0.0) {
        m_r0 = gammalib::rad2deg;
    }

    // Precompute
    m_w[0] = 1.0 / m_r0;
    m_w[1] = m_pv[1]*m_pv[1] + m_pv[2]*m_pv[2];
    m_w[2] = m_w[1] + 1.0;
    m_w[3] = m_w[1] - 1.0;
    
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
 * Deproject pixel (x,y) coordinates in the plane of projection to native
 * spherical coordinates (phi,theta).
 *
 * This method has been adapted from the wcslib function prj.c::sinx2s().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previously (as this will be done by the constructor).
 ***************************************************************************/
void GWcsSIN::prj_x2s(int nx, int ny, int sxy, int spt, 
                      const double* x, const double* y,
                      double* phi, double* theta, int* stat) const
{
    // Set tolerance
    const double tol = 1.0e-13;

    // Initialize projection if required
    if (!m_prjset) {
        prj_set();
    }

    // Initialise status code and statistics
    int status    = 0;
    int n_invalid = 0;

    // Get obliqueness parameters
    double xi  = m_pv[1];
    double eta = m_pv[2];

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
        double  x0   = (*xp + m_x0) * m_w[0];
        double* phip = phi + rowoff;
        for (int iy = 0; iy < my; ++iy, phip += rowlen) {
            *phip = x0;
        }
    }

    // Do y dependence
    const double* yp     = y;
    double*       phip   = phi;
    double*       thetap = theta;
    int*          statp  = stat;
    for (int iy = 0; iy < ny; ++iy, yp += sxy) {
        double y0  = (*yp + m_y0) * m_w[0];
        double y02 = y0*y0;
        for (int ix = 0; ix < mx; ++ix, phip += spt, thetap += spt) {

            // Compute intermediaries
            double x0 = *phip;
            double r2 = x0*x0 + y02;

            // Orthographic projection
            if (m_w[1] == 0.0) {

                // Compute Phi
                if (r2 != 0.0) {
                    *phip = gammalib::atan2d(x0, -y0);
                }
                else {
                    *phip = 0.0;
                }

                // Compute Theta
                if (r2 < 0.5) {
                    *thetap = gammalib::acosd(std::sqrt(r2));
                }
                else if (r2 <= 1.0) {
                    *thetap = gammalib::asind(std::sqrt(1.0 - r2));
                }
                else {
                    *thetap    = 0.0;
                    *(statp++) = 1;
                    status     = 3;
                    n_invalid++;
                    continue;
                }

            } // endif: orthographic projection
            
            // "Synthesis" projection
            else {

                // ...
                double z;
                double xy = x0*xi + y0*eta;

                // Use small angle formula
                if (r2 < 1.0e-10) {
                    z       = r2/2.0;
                    *thetap = 90.0 - gammalib::rad2deg*std::sqrt(r2/(1.0 + xy));

                }
                else {
                    double a = m_w[2];
                    double b = xy - m_w[1];
                    double c = r2 - xy - xy + m_w[3];
                    double d = b*b - a*c;

                    // Check for a solution
                    if (d < 0.0) {
                        *phip      = 0.0;
                        *thetap    = 0.0;
                        *(statp++) = 1;
                        status     = 3;
                        n_invalid++;
                        continue;
                    }

                    // ...
                    d = std::sqrt(d);

                    // Choose solution closest to pole
                    double sinth1 = (-b + d)/a;
                    double sinth2 = (-b - d)/a;
                    double sinthe = (sinth1 > sinth2) ? sinth1 : sinth2;
                    if (sinthe > 1.0) {
                        if (sinthe-1.0 < tol) {
                            sinthe = 1.0;
                        }
                        else {
                            sinthe = (sinth1 < sinth2) ? sinth1 : sinth2;
                        }
                    }

                    // ...
                    if (sinthe < -1.0) {
                        if (sinthe+1.0 > -tol) {
                            sinthe = -1.0;
                        }
                    }

                    // ...
                    if (sinthe > 1.0 || sinthe < -1.0) {
                        *phip      = 0.0;
                        *thetap    = 0.0;
                        *(statp++) = 1;
                        status     = 3;
                        n_invalid++;
                        continue;
                    }

                    // ...
                    *thetap = gammalib::asind(sinthe);
                    z = 1.0 - sinthe;
                }

                // ...
                double x1 = -y0 + eta * z;
                double y1 =  x0 -  xi * z;
                if (x1 == 0.0 && y1 == 0.0) {
                    *phip = 0.0;
                }
                else {
                    *phip = gammalib::atan2d(y1,x1);
                }
                
            } // endelse: Synthesis projection

            // Reset error
            *(statp++) = 0;
            
        } // endfor
        
    } // endfor: do Y dependence

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
 * Project native spherical coordinates (phi,theta) to pixel (x,y)
 * coordinates in the plane of projection.
 *
 * This method has been adapted from the wcslib function prj.c::sins2x().
 * The interface follows very closely that of wcslib. In contrast to the
 * wcslib routine, however, the method assumes that the projection has been
 * setup previously (as this will be done by the constructor).
 ***************************************************************************/
void GWcsSIN::prj_s2x(int nphi, int ntheta, int spt, int sxy,
                      const double* phi, const double* theta,
                      double* x, double* y, int* stat) const
{
    // Initialize projection if required
    if (!m_prjset) {
        prj_set();
    }

    // Initialise status code and statistics
    int status    = 0;
    int n_invalid = 0;

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

        // ...
        double z;
        double costhe;
        double t = (90.0 - std::abs(*thetap)) * gammalib::deg2rad;
        if (t < 1.0e-5) {
            if (*thetap > 0.0) {
                z = t * t / 2.0;
            }
            else {
                z = 2.0 - t * t / 2.0;
            }
            costhe = t;
        }
        else {
            z      = 1.0 - gammalib::sind(*thetap);
            costhe = gammalib::cosd(*thetap);
        }
        double r = m_r0 * costhe;

        // Orthographic projection
        if (m_w[1] == 0.0) {
            int istat = 0;
            if (m_bounds) {
                if (*thetap < 0.0) {
                    istat  = 1;
                    status = 4;
                    n_invalid++;
                }
            }
            for (int iphi = 0; iphi < mphi; ++iphi, xp += sxy, yp += sxy) {
                *xp =  r*(*xp) - m_x0;
                *yp = -r*(*yp) - m_y0;
                *(statp++) = istat;
            }
        }

        // Synthesis projection
        else {
        
            // ...
            z *= m_r0;
            double z1 = m_pv[1] * z - m_x0;
            double z2 = m_pv[2] * z - m_y0;

            // ...
            for (int iphi = 0; iphi < mphi; ++iphi, xp += sxy, yp += sxy) {
                int istat = 0;
                if (m_bounds) {
                    double t = -gammalib::atand(m_pv[1]*(*xp) - m_pv[2]*(*yp));
                    if (*thetap < t) {
                        istat  = 1;
                        status = 4;
                        n_invalid++;
                    }
                }
                *xp =  r * (*xp) + z1;
                *yp = -r * (*yp) + z2;
                *(statp++) = istat;
            }
            
        } // endelse: Synthesis projection
        
    } // endfor: do theta dependence

    // Handle status code
    if (status == 4) {
        throw GException::wcs_invalid_phi_theta(G_PRJ_S2X, n_invalid);
    }
    
    // Return
    return;
}
