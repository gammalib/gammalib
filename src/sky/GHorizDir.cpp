/***************************************************************************
 *                  GHorizDir.cpp - Horizontal direction class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @brief Horiz direction class implementation
 * @author K. Kosack
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
 * @param[in] dir Horiz direction.
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
 * @param[in] dir Horiz direction.
 * @return Horiz direction.
 ***************************************************************************/
GHorizDir& GHorizDir::operator= (const GHorizDir& dir)
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
 * @brief Clear sky direction
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
 * @brief Clone sky direction
 *
 * @return Pointer to deep copy of sky direction.
 ***************************************************************************/
GHorizDir* GHorizDir::clone(void) const
{
    // Clone sky direction
    return new GHorizDir(*this);
}


/***********************************************************************//**
 * @brief Set equatorial sky direction (radians)
 *
 * @param[in] ra Right Ascension in radians.
 * @param[in] dec Declination in radians.
 *
 * Sets Right Ascension and Declination in radians.
 ***************************************************************************/
void GHorizDir::altaz(const double& alt, const double& az)
{

    // Set direction
    m_alt  = alt;
    m_az = az;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set equatorial sky direction (degrees)
 *
 * @param[in] ra Right Ascension in degrees.
 * @param[in] dec Declination in degrees.
 *
 * Sets Right Ascension and Declination in degrees.
 ***************************************************************************/
void GHorizDir::altaz_deg(const double& alt, const double& az)
{

    // Set direction
    m_ra  = alt  * gammalib::deg2rad;
    m_dec = az * gammalib::deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky direction from 3D vector in horizontal coordinates
 *
 * @param[in] vector 3D vector.
 *
 * Convert a 3-dimensional vector in celestial coordinates into a sky
 * direction. The transformation is given by
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


    // Convert vector into sky position
    m_alt = std::asin(vector[2]);
    m_az  = std::atan2(vector[1], vector[0]);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Rotate sky direction by zenith and azimuth angle
 *
 * @param[in] phi Azimuth angle (deg).
 * @param[in] theta Zenith angle (deg).
 *
 * Rotate sky direction by a zenith and azimuth angle given in the system
 * of the sky direction and aligned in celestial coordinates. 
 * The azimuth angle is counted counter clockwise from celestial north
 * (this is identical to the astronomical definition of a position angle).
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

    // Rotate vector into celestial coordinates
    GVector dir = rot * native;

    // Convert vector into sky position
    celvector(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return sky direction as 3D vector in celestial coordinates
 *
 * @return Horiz direction as 3D vector in celestial coordinates.
 ***************************************************************************/
GVector GHorizDir::celvector(void) const
{
    // Compute 3D vector
    double  cosaz  = std::cos(m_az);
    double  sinaz  = std::sin(m_az);
    double  cosalt = std::cos(m_alt);
    double  sinalt = std::sin(m_alt);
    GVector vector(cosdec*cosaz, cosdec*sinaz, sinalt);

    // Return vector
    return vector;
}


/***********************************************************************//**
 * @brief Compute angular distance between sky directions in radians
 *
 * @param[in] dir Horiz direction to which distance is to be computed.
 * @return Angular distance in radians.
 *
 * Computes the angular distance between two sky directions in radians.
 ***************************************************************************/
double GHorizDir::dist(const GHorizDir& dir) const
{
    // Initialise cosine of distance
    double cosdis;

    // Compute dependent on coordinate system availability. This speeds
    // up things by avoiding unnecessary coordinate transformations.

    cosdis = sin(m_alt) * sin(dir.alt()) +
      cos(m_alt) * cos(dir.alt()) *
      std::cos(dir.az() - m_az);
    
    // Compute distance (use argument save GTools function)
    double dist = gammalib::acos(cosdis);

    // Return distance
    return dist;
}


/***********************************************************************//**
 * @brief Compute angular distance between sky directions in degrees
 *
 * @param[in] dir Horiz direction to which distance is to be computed.
 * @return Angular distance in degrees.
 *
 * Computes the angular distance between two sky directions in degrees.
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
 * @return String containing sky direction information.
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
    m_az        = 0.0;
    m_alt       = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Horiz direction.
 ***************************************************************************/
void GHorizDir::copy_members(const GHorizDir& dir)
{
    m_az        = dir.m_az;
    m_alt       = dir.m_alt;

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

/***********************************************************************
 * @brief Equality operator
 *
 * @param[in] a First horiz direction.
 * @param[in] b Second horiz direction.
 *
 * Compare two horiz directions. If the coordinate is at the pole, the
 * Azimuth value is irrelevant.
 *
 * @todo: really should test for being within some tolerance here
 *
 ***************************************************************************/
bool operator==(const GHorizDir &a, const GHorizDir &b)
{
    // Initialise result
    bool equal = false;
    
    // Compute dependent on coordinate system availability. This speeds
    // up things by avoiding unnecessary coordinate transformations.
    
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
 * @param[in] a First sky direction.
 * @param[in] b Second sky direction.
 ***************************************************************************/
bool operator!=(const GHorizDir &a, const GHorizDir &b)
{
    // Return non equality
    return (!(a==b));
}

