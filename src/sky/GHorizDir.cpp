/***************************************************************************
 *                 GHorizDir.cpp - Horizontal direction class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Karl Kosack                                      *
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
 * @file GHorizDir.cpp
 * @brief Horizontal direction class implementation
 * @author Karl Kosack
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GHorizDir.hpp"
#include "GMatrix.hpp"
#include "GVector.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                          Constructors/destructors                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GHorizDir::GHorizDir(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Horizontal direction.
 ***************************************************************************/
GHorizDir::GHorizDir(const GHorizDir& dir)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GHorizDir::~GHorizDir(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] dir Horizontal direction.
 * @return Horizontal direction.
 ***************************************************************************/
GHorizDir& GHorizDir::operator=(const GHorizDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(dir);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear horizontal direction
 ***************************************************************************/
void GHorizDir::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone horizontal direction
 *
 * @return Pointer to deep copy of horizontal direction.
 ***************************************************************************/
GHorizDir* GHorizDir::clone(void) const
{
    // Clone horizon direction
    return new GHorizDir(*this);
}


/***********************************************************************//**
 * @brief Set horizontal direction (radians)
 *
 * @param[in] alt Altitude in radians.
 * @param[in] az Azimuth in radians.
 *
 * Sets Altitude and Azimuth in radians.
 ***************************************************************************/
void GHorizDir::altaz(const double& alt, const double& az)
{

    // Set direction
    m_alt = alt;
    m_az  = az;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set horizontal direction (degrees)
 *
 * @param[in] alt Altitude in radians.
 * @param[in] az Azimuth in radians.
 *
 * Sets Altitude and Azimuth in degrees.
 ***************************************************************************/
void GHorizDir::altaz_deg(const double& alt, const double& az)
{

    // Set direction
    m_az  = az  * gammalib::deg2rad;
    m_alt = alt * gammalib::deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set horizontal direction from 3D vector
 *
 * @param[in] vector 3D vector.
 *
 * Convert a 3-dimensional vector in into a horizontal direction. The
 * transformation is given by
 *
 * \f[
 *    \alpha = \arctan \left( \frac{x_1}{x_0} \right)
 * \f]
 *
 * \f[
 *    \delta = \arcsin x_2
 * \f]
 ***************************************************************************/
void GHorizDir::celvector(const GVector& vector)
{
    // Convert vector into horizontal direction
    m_alt = std::asin(vector[2]);
    m_az  = std::atan2(vector[1], vector[0]);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Rotate horizontal direction by zenith and azimuth angle
 *
 * @param[in] phi Azimuth angle (deg).
 * @param[in] theta Zenith angle (deg).
 *
 * Rotate horizontal direction by a zenith and azimuth angle given in the
 * system of the horizontal direction and aligned in TBD. 
 * The azimuth angle is counted counter clockwise from TBD.
 ***************************************************************************/
void GHorizDir::rotate_deg(const double& phi, const double& theta)
{

    // Allocate Euler and rotation matrices
    GMatrix ry;
    GMatrix rz;

    // Set up rotation matrix to rotate from native coordinates to
    // celestial coordinates
    ry.eulery(m_alt * gammalib::rad2deg - 90.0);
    rz.eulerz(-m_az * gammalib::rad2deg);
    GMatrix rot = (ry * rz).transpose();

    // Set up native coordinate vector
    double phi_rad   = phi   * gammalib::deg2rad;
    double theta_rad = theta * gammalib::deg2rad;
    double cos_phi   = std::cos(phi_rad);
    double sin_phi   = std::sin(phi_rad);
    double cos_theta = std::cos(theta_rad);
    double sin_theta = std::sin(theta_rad);
    GVector native(-cos_phi*sin_theta, sin_phi*sin_theta, cos_theta);

    // Rotate vector into horizontal coordinates
    GVector dir = rot * native;

    // Convert vector into horizontal direction
    celvector(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return horizontal direction as 3D vector
 *
 * @return Horizontal direction as 3D vector.
 ***************************************************************************/
GVector GHorizDir::celvector(void) const
{
    // Compute 3D vector
    double  cosaz  = std::cos(m_az);
    double  sinaz  = std::sin(m_az);
    double  cosalt = std::cos(m_alt);
    double  sinalt = std::sin(m_alt);
    GVector vector(cosalt*cosaz, cosalt*sinaz, sinalt);

    // Return vector
    return vector;
}


/***********************************************************************//**
 * @brief Compute angular distance to horizontal direction in radians
 *
 * @param[in] dir Horizontal direction.
 * @return Angular distance in radians.
 *
 * Computes the angular distance to a specified horizontal direction in
 * radians.
 ***************************************************************************/
double GHorizDir::dist(const GHorizDir& dir) const
{
    // Compute distance
    double cosdis = std::sin(m_alt) * std::sin(dir.alt()) +
                    std::cos(m_alt) * std::cos(dir.alt()) *
                    std::cos(dir.az() - m_az);
    
    // Compute distance (use argument save GTools function)
    double dist = gammalib::acos(cosdis);

    // Return distance
    return dist;
}


/***********************************************************************//**
 * @brief Compute angular distance to horizontal direction in degrees
 *
 * @param[in] dir Horizontal direction.
 * @return Angular distance in degrees.
 *
 * Computes the angular distance to a specified horizontal direction in
 * degrees.
 ***************************************************************************/
double GHorizDir::dist_deg(const GHorizDir& dir) const
{
    // Return distance in degrees
    return (dist(dir) * gammalib::rad2deg);
}



/***********************************************************************//**
 * @brief Print horizontal direction information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing horizontal direction information.
 ***************************************************************************/
std::string GHorizDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

      // Put coordinates in string
      result = "(Az,Alt)=("+gammalib::str(m_az*gammalib::rad2deg) + "," +
        gammalib::str(m_alt*gammalib::rad2deg)+")";
      
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
void GHorizDir::init_members(void)
{
    // Initialise members
    m_az  = 0.0;
    m_alt = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Horizontal direction.
 ***************************************************************************/
void GHorizDir::copy_members(const GHorizDir& dir)
{
    // Copy members
    m_az  = dir.m_az;
    m_alt = dir.m_alt;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GHorizDir::free_members(void)
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
 * @param[in] a First horizontal direction.
 * @param[in] b Second horizontal direction.
 *
 * Compare two horizontal directions. If the coordinate is at the pole,
 * the azimuth value is irrelevant.
 *
 * @todo: really should test for being within some tolerance here
 ***************************************************************************/
bool operator==(const GHorizDir &a, const GHorizDir &b)
{
    // Initialise result
    bool equal = false;
    
    // Do check
    if (std::abs(std::abs(a.m_alt) - 90.0) < 1e-10) {
        equal = (a.m_alt == b.m_alt);
    }
    else {
        equal = (a.m_alt == b.m_alt && a.m_az == b.m_az);
    }
    
    // Return equality
    return equal;
}


/***********************************************************************//**
 * @brief Non equality operator
 *
 * @param[in] a First horizontal direction.
 * @param[in] b Second horizontal direction.
 ***************************************************************************/
bool operator!=(const GHorizDir &a, const GHorizDir &b)
{
    // Return non equality
    return (!(a==b));
}
