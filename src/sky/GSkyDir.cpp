/***************************************************************************
 *                     GSkyDir.cpp - Sky direction class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2019 by Juergen Knoedlseder                         *
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
 * @file GSkyDir.cpp
 * @brief Sky direction class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GSkyDir.hpp"
#include "GTime.hpp"
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
GSkyDir::GSkyDir(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Sky direction.
 ***************************************************************************/
GSkyDir::GSkyDir(const GSkyDir& dir)
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
GSkyDir::~GSkyDir(void)
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
 * @param[in] dir Sky direction.
 * @return Sky direction.
 ***************************************************************************/
GSkyDir& GSkyDir::operator= (const GSkyDir& dir)
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
void GSkyDir::clear(void)
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
GSkyDir* GSkyDir::clone(void) const
{
    // Clone sky direction
    return new GSkyDir(*this);
}


/***********************************************************************//**
 * @brief Set equatorial sky direction (radians)
 *
 * @param[in] ra Right Ascension in radians.
 * @param[in] dec Declination in radians.
 *
 * Sets Right Ascension and Declination in radians.
 ***************************************************************************/
void GSkyDir::radec(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = false;
    m_has_radec = true;
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = false;
    m_has_radec_cache = false;
    #endif

    // Set direction
    m_ra  = ra;
    m_dec = dec;

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
void GSkyDir::radec_deg(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = false;
    m_has_radec = true;
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = false;
    m_has_radec_cache = false;
    #endif

    // Set direction
    m_ra  = ra  * gammalib::deg2rad;
    m_dec = dec * gammalib::deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set galactic sky direction (radians)
 *
 * @param[in] l Galactic longitude in radians.
 * @param[in] b Galactic latitude in radians.
 *
 * Sets Galactic longitude and latitude in radians.
 ***************************************************************************/
void GSkyDir::lb(const double& l, const double& b)
{
    // Set attributes
    m_has_lb    = true;
    m_has_radec = false;
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = false;
    m_has_radec_cache = false;
    #endif

    // Set direction
    m_l = l;
    m_b = b;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set galactic sky direction (degrees)
 *
 * @param[in] l Galactic longitude in degrees.
 * @param[in] b Galactic latitude in degrees.
 *
 * Sets Galactic longitude and latitude in degrees.
 ***************************************************************************/
void GSkyDir::lb_deg(const double& l, const double& b)
{
    // Set attributes
    m_has_lb    = true;
    m_has_radec = false;
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = false;
    m_has_radec_cache = false;
    #endif

    // Set direction
    m_l = l * gammalib::deg2rad;
    m_b = b * gammalib::deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky direction from 3D vector in celestial coordinates
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
void GSkyDir::celvector(const GVector& vector)
{
    // Set attributes
    m_has_lb    = false;
    m_has_radec = true;
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = false;
    m_has_radec_cache = false;
    #endif

    // Convert vector into sky position
    m_dec = std::asin(vector[2]);
    m_ra  = std::atan2(vector[1], vector[0]);

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
void GSkyDir::rotate_deg(const double& phi, const double& theta)
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb) {
        gal2equ();
    }

    // Allocate Euler and rotation matrices
    GMatrix ry;
    GMatrix rz;

    // Set up rotation matrix to rotate from native coordinates to
    // celestial coordinates
    ry.eulery(m_dec * gammalib::rad2deg - 90.0);
    rz.eulerz(-m_ra * gammalib::rad2deg);
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
 * @brief Precess sky direction
 *
 * @param[in] from_epoch Epoch of the current coordinate.
 * @param[in] to_epoch Epoch of the returned precessed coordinate.
 *
 * Precesses the sky direction from one epoch to another.
 *
 * The method is adapted from a set of fortran subroutines based on
 * (a) pages 30-34 of the Explanatory Supplement to the AE,
 * (b) Lieske, et al. (1977) A&A 58, 1-16, and
 * (c) Lieske (1979) A&A 73, 282-284.
 ***************************************************************************/
void GSkyDir::precess(const double& from_epoch, const double& to_epoch)
{
    // Set constants
    const double arcsec = gammalib::pi / 180.0 / 3600.0;

    // Continue only if epochs differ
    if (from_epoch != to_epoch) {

        // t0, t below correspond to Lieske's big T and little T
        double t0  = (from_epoch - 2000.0)     / 100.0;
        double t   = (to_epoch   - from_epoch) / 100.0;
        double t02 = t0*t0;
        double t2  = t*t;
        double t3  = t2*t;

        // a,b,c below correspond to Lieske's zeta_A, z_A and theta_A
        double a = ((2306.2181 + 1.39656  * t0 - 0.000139 * t02) * t +
                    (0.30188   - 0.000344 * t0) * t2 + 0.017998 * t3) * arcsec;
        double b = ((2306.2181 + 1.39656  * t0 - 0.000139 * t02) * t +
                    (1.09468   + 0.000066 * t0) * t2 + 0.018203 * t3) * arcsec;
        double c = ((2004.3109 - 0.85330  * t0 - 0.000217 * t02) * t +
                    (-0.42665  - 0.000217 * t0) * t2 - 0.041833 * t3) * arcsec;

        // Compute sines and cosines
        double sina = std::sin(a);
        double cosa = std::cos(a);
        double sinb = std::sin(b);
        double cosb = std::cos(b);
        double sinc = std::sin(c);
        double cosc = std::cos(c);

        // Setup precession rotation matrix
        double xx =  cosa*cosc*cosb - sina*sinb;
        double yx = -sina*cosc*cosb - cosa*sinb;
        double zx = -sinc*cosb;
        double xy =  cosa*cosc*sinb + sina*cosb;
        double yy = -sina*cosc*sinb + cosa*cosb;
        double zy = -sinc*sinb;
        double xz =  cosa*sinc;
        double yz = -sina*sinc;
        double zz =  cosc;

        // Get vector of sky position
        GVector vector = celvector();

        // Perform the rotation:
        double x2 = xx*vector[0] + yx*vector[1] + zx*vector[2];
        double y2 = xy*vector[0] + yy*vector[1] + zy*vector[2];
        double z2 = xz*vector[0] + yz*vector[1] + zz*vector[2];

        // Setup rotated vector
        GVector rotated_vector(x2, y2, z2);

        // Transform vector to sky direction
        celvector(rotated_vector);

    } // endif: epochs differed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky direction to direction of Sun
 *
 * @param[in] time Time.
 * @param[in] epoch Julian epoch.
 *
 * Sets the sky direction to the direction of the Sun at a given @p time.
 *
 * The equations were taken from the Astronomical Almanac. They calculate the
 * apparent coordinates of the Sun to a precision of about 0.01 degrees
 * (36 arcsec), for dates between 1950 and 2050.
 *
 * See https://en.wikipedia.org/wiki/Position_of_the_Sun
 ***************************************************************************/
void GSkyDir::sun(const GTime& time, const double& epoch)
{
    // Compute number of days since Greenwich noon, Terrestrial Time, on
    // 1 January 2000
    double n = time.jd() - 2451545.0;

    // Compute obliquity of the ecliptic in radians
    double eps_rad = (23.439 - 0.0000004 * n) * gammalib::deg2rad;

    // Compute mean longitude of the Sun in degrees, corrected for the
    // aberration of light. Put the mean longitude in the interval [0,360[
    double mean_longitude = 280.460 + 0.9856474 * n;
    while (mean_longitude < 0.0) {
        mean_longitude += 360.0;
    }
    while (mean_longitude >= 360.0) {
        mean_longitude -= 360.0;
    }
    
    // Compute the mean anomaly of the Sun in radians
    double mean_anomaly = (357.528 + 0.9856003 * n) * gammalib::deg2rad;

    // Compute the ecliptic longitude of the Sun in degrees
    double ecliptic_longitude = mean_longitude +
                                1.915 * std::sin(mean_anomaly) +
                                0.020 * std::sin(2.0*mean_anomaly);

    // Compute sine and cosine
    double ecliptic_longitude_rad = ecliptic_longitude * gammalib::deg2rad;
    double sin_ecliptic_longitude = std::sin(ecliptic_longitude_rad);
    double cos_ecliptic_longitude = std::cos(ecliptic_longitude_rad);

    // Compute Right Ascension and Declination of the Sun in degrees
    double ra  = std::atan2(std::cos(eps_rad) * sin_ecliptic_longitude,
                            cos_ecliptic_longitude) * gammalib::rad2deg;
    double dec = std::asin(std::sin(eps_rad) * sin_ecliptic_longitude) *
                           gammalib::rad2deg;

    // Set position
    radec_deg(ra, dec);

    // Precess sky position to requested epoch
    precess(time.julian_epoch(), epoch);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky direction to direction of Moon
 *
 * @param[in] time Time.
 * @param[in] epoch Julian epoch.
 *
 * Sets the sky direction to the direction of the Moon at a given @p time.
 *
 * The equations are derived from the Chapront ELP2000/82 Lunar Theory
 * (Chapront-Touze' and Chapront, 1983, 124, 50), as described by Jean Meeus
 * in Chapter 47 of ``Astronomical Algorithms'' (Willmann-Bell, Richmond),
 * 2nd edition, 1998. Meeus quotes an approximate accuracy of 10 arcsec in
 * longitude and 4 arcsec in latitude, but he does not give the time range
 * for this accuracy.
 *
 * The method was compared to the results obtained using PyEphem for the
 * period 2000 - 2050. Differences in Right Ascension and Declination are
 * always smaller than 0.015 degrees, corresponding to about 1 arcmin.
 * The differences increase over time, meaning that the largest distances
 * are reached by the end of the period.
 ***************************************************************************/
void GSkyDir::moon(const GTime& time, const double& epoch)
{
    // Set longitude constants. Note that the cos_lng contants are commented
    // out since they are actually not used. They could be useful in case
    // that the distance to the Moon should be computed, therefore we leave
    // them in comments in the code.
    const int d_lng[]   = { 0, 2, 2, 0, 0, 0, 2, 2, 2, 2,
                            0, 1, 0, 2, 0, 0, 4, 0, 4, 2,
                            2, 1, 1, 2, 2, 4, 2, 0, 2, 2,
                            1, 2, 0, 0, 2, 2, 2, 4, 0, 3,
                            2, 4, 0, 2, 2, 2, 4, 0, 4, 1,
                            2, 0, 1, 3, 4, 2, 0, 1, 2, 2};
    const int m_lng[]   = { 0, 0, 0, 0, 1, 0, 0,-1, 0,-1,
                            1, 0, 1, 0, 0, 0, 0, 0, 0, 1,
                            1, 0, 1,-1, 0, 0, 0, 1, 0,-1,
                            0,-2, 1, 2,-2, 0, 0,-1, 0, 0,
                            1,-1, 2, 2, 1,-1, 0, 0,-1, 0,
                            1, 0, 1, 0, 0,-1, 2, 1, 0,0};
    const int mp_lng[]  = { 1,-1, 0, 2, 0, 0,-2,-1, 1, 0,
                           -1, 0, 1, 0, 1, 1,-1, 3,-2,-1,
                            0,-1, 0, 1, 2, 0,-3,-2,-1,-2,
                            1, 0, 2, 0,-1, 1, 0,-1, 2,-1,
                            1,-2,-1,-1,-2, 0, 1, 4, 0,-2,
                            0, 2, 1,-2,-3, 2, 1,-1, 3,-1};
    const int f_lng[]   = { 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
                            0, 0, 0,-2, 2,-2, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 2, 0,
                            0, 0, 0, 0, 0,-2, 2, 0, 2, 0,
                            0, 0, 0, 0, 0,-2, 0, 0, 0, 0,
                           -2,-2, 0, 0, 0, 0, 0, 0, 0,-2};
    const int sin_lng[] = {  6288774, 1274027,  658314, 213618,-185116,
                             -114332,   58793,   57066,  53322,  45758,
                              -40923,  -34720,  -30383,  15327, -12528,
                               10980,   10675,   10034,   8548,  -7888,
                               -6766,   -5163,    4987,   4036,   3994,
                                3861,    3665,   -2689,  -2602,   2390,
                               -2348,    2236,   -2120,  -2069,   2048,
                               -1773,   -1595,    1215,  -1110,   -892,
                                -810,     759,    -713,   -700,    691,
                                 596,     549,     537,    520,   -487,
                                -399,    -381,     351,   -340,    330,
                                 327,    -323,     299,    294,      0};
    /*
    const int cos_lng[] = {-20905355,-3699111,-2955968,-569925,  48888,
                               -3149,  246158, -152138,-170733,-204586,
                             -129620,  108743,  104755,  10321,      0,
                               79661,  -34782,  -23210, -21636,  24208,
                               30824,   -8379,  -16675, -12831, -10445,
                              -11650,   14403,   -7003,      0,  10056,
                                6322,   -9884,    5751,      0,  -4950,
                                4130,       0,   -3958,      0,   3258,
                                2616,   -1897,   -2117,   2354,      0,
                                   0,   -1423,   -1117,  -1571,  -1739,
                                   0,   -4421,       0,      0,      0,
                                   0,    1165,       0,      0,   8752};
    */

    // Set latitude constants
    const int d_lat[]   = { 0, 0, 0, 2, 2, 2, 2, 0, 2, 0,
                            2, 2, 2, 2, 2, 2, 2, 0, 4, 0,
                            0, 0, 1, 0, 0, 0, 1, 0, 4, 4,
                            0, 4, 2, 2, 2, 2, 0, 2, 2, 2,
                            2, 4, 2, 2, 0, 2, 1, 1, 0, 2,
                            1, 2, 0, 4, 4, 1, 4, 1, 4, 2};
    const int m_lat[]   = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           -1, 0, 0, 1,-1,-1,-1, 1, 0, 1,
                            0, 1, 0, 1, 1, 1, 0, 0, 0, 0,
                            0, 0, 0, 0,-1, 0, 0, 0, 0, 1,
                            1, 0,-1,-2, 0, 1, 1, 1, 1, 1,
                            0,-1, 1, 0,-1, 0, 0, 0,-1,-2};
    const int mp_lat[]  = { 0, 1, 1, 0,-1,-1, 0, 2, 1, 2,
                            0,-2, 1, 0,-1, 0,-1,-1,-1, 0,
                            0,-1, 0, 1, 1, 0, 0, 3, 0,-1,
                            1,-2, 0, 2, 1,-2, 3, 2,-3,-1,
                            0, 0, 1, 0, 1, 1, 0, 0,-2,-1,
                            1,-2, 2,-2,-1, 1, 1,-1, 0, 0};
    const int f_lat[]   = { 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,
                           -1,-1, 1,-1, 1, 1,-1,-1,-1, 1,
                            3, 1, 1, 1,-1,-1,-1, 1,-1, 1,
                           -3, 1,-3,-1,-1, 1,-1, 1,-1, 1,
                            1, 1, 1,-1, 3,-1,-1, 1,-1,-1,
                            1,-1, 1,-1,-1,-1,-1,-1,-1, 1};
    const int sin_lat[] = {5128122, 280602, 277693, 173237, 55413,
                             46271,  32573,  17198,   9266,  8822,
                              8216,   4324,   4200,  -3359,  2463,
                              2211,   2065,  -1870,   1828, -1794,
                             -1749,  -1565,  -1491,  -1475, -1410,
                             -1344,  -1335,   1107,   1021,   833,
                               777,    671,    607,    596,   491,
                              -451,    439,    422,    421,  -366,
                              -351,    331,    315,    302,  -283,
                              -229,    223,    223,   -220,  -220,
                              -185,    181,   -177,    176,   166,
                              -164,    132,   -119,    115,   107};

    // Set constants for computation of nutation. The constants have been
    // taken from https://idlastro.gsfc.nasa.gov/ftp/pro/astro/nutate.pro
    // which is based on Chapter 22 of "Astronomical Algorithms" by Jean
    // Meeus (1998, 2nd ed.) which is based on the 1980 IAU Theory of
    // Nutation and includes all terms larger than 0.0003 arcsec.
    const int nutate_d_lng[]    = { 0,-2, 0, 0, 0, 0,-2, 0, 0,-2,
                                   -2,-2, 0, 2, 0, 2, 0, 0,-2, 0,
                                    2, 0, 0,-2, 0,-2, 0, 0, 2,-2,
                                    0,-2, 0, 0, 2, 2, 0,-2, 0, 2,
                                    2,-2,-2, 2, 2, 0,-2,-2, 0,-2,
                                   -2, 0,-1,-2, 1, 0, 0,-1, 0, 0,
                                    2, 0, 2};
    const int nutate_m_lng[]    = { 0, 0, 0, 0, 1, 0, 1, 0, 0,-1,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 2, 0, 2,
                                    1, 0,-1, 0, 0, 0, 1, 1,-1, 0,
                                    0, 0, 0, 0, 0,-1,-1, 0, 0, 0,
                                    1, 0, 0, 1, 0, 0, 0,-1, 1,-1,
                                   -1, 0,-1};
    const int nutate_mp_lng[]   = { 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
                                    1, 0,-1, 0, 1,-1,-1, 1, 2,-2,
                                    0, 2, 2, 1, 0, 0,-1, 0,-1, 0,
                                    0, 1, 0, 2,-1, 1, 0, 1, 0, 0,
                                    1, 2, 1,-2, 0, 1, 0, 0, 2, 2,
                                    0, 1, 1, 0, 0, 1,-2, 1, 1, 1,
                                   -1, 3, 0};
    const int nutate_f_lng[]    = { 0, 2, 2, 0, 0, 0, 2, 2, 2, 2,
                                    0, 2, 2, 0, 0, 2, 0, 2, 0, 2,
                                    2, 2, 0, 2, 2, 2, 2, 0, 0, 2,
                                    0, 0, 0,-2, 2, 2, 2, 0, 2, 2,
                                    0, 2, 2, 0, 0, 0, 2, 0, 2, 0,
                                    2,-2, 0, 0, 0, 2, 2, 0, 0, 2,
                                    2, 2, 2};
    const int nutate_om_lng[]   = { 1, 2, 2, 2, 0, 0, 2, 1, 2, 2,
                                    0, 1, 2, 0, 1, 2, 1, 1, 0, 1,
                                    2, 2, 0, 2, 0, 0, 1, 0, 1, 2,
                                    1, 1, 1, 0, 1, 2, 2, 0, 2, 1,
                                    0, 2, 1, 1, 1, 0, 1, 1, 1, 1,
                                    1, 0, 0, 0, 0, 0, 2, 0, 0, 2,
                                    2, 2, 2};
    const int nutate_sin_lng[] = {-171996, -13187, -2274, 2062, 1426,
                                      712,   -517,  -386, -301,  217,
                                     -158,    129,   123,   63,   63,
                                      -59,    -58,   -51,   48,   46,
                                      -38,    -31,    29,   29,   26,
                                      -22,     21,    17,   16,  -16,
                                      -15,    -13,   -12,   11,  -10,
                                       -8,      7,    -7,   -7,   -7,
                                        6,      6,     6,   -6,   -6,
                                        5,     -5,    -5,   -5,    4,
                                        4,      4,    -4,   -4,   -4,
                                        3,     -3,    -3,   -3,   -3,
                                       -3,     -3,    -3};
    const int nutate_cos_lng[] = {  92025,   5736,   977, -895,   54,
                                       -7,    224,   200,  129,  -95,
                                        0,    -70,   -53,    0,  -33,
                                       26,     32,    27,    0,  -24,
                                       16,     13,     0,  -12,    0,
                                        0,    -10,     0,   -8,    7,
                                        9,      7,     6,    0,    5,
                                        3,     -3,     0,    3,    3,
                                        0,     -3,    -3,    3,    3,
                                        0,      3,     3,    3,    0,
                                        0,      0,     0,    0,    0,
                                        0,      0,     0,    0,    0,
                                        0,      0,     0};
    const double nutate_sdelt[] = {-174.2, -1.6, -0.2,  0.2, -3.4,
                                      0.1,  1.2, -0.4,  0.0, -0.5,
                                      0.0,  0.1,  0.0,  0.0,  0.1,
                                      0.0, -0.1,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0, -0.1,  0.0,  0.1,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0};
    const double nutate_cdelt[] =    {8.9, -3.1, -0.5,  0.5, -0.1,
                                      0.0, -0.6,  0.0, -0.1,  0.3,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0,  0.0,  0.0,
                                      0.0,  0.0,  0.0};

    // Compute Julian centuries from 1900.0
    double t = (time.jd() - 2451545.0) / 36525.0;

    // Compute Moon's mean longitude referred to mean equinox of the date
    // in degrees and radians
    double Lprimed = 218.3164477 +
                     (481267.88123421 +
                     (-0.0015786 +
                     (1.0/538841.0 - 1.0/65194000.0 * t) * t) * t) * t;
    double Lprime = Lprimed * gammalib::deg2rad;

    // Compute Moon's mean elongation in degrees
    double D = (297.8501921 +
               (445267.1114034 +
               (-0.0018819 +
               (1.0/545868.0 - 1.0/113065000.0 * t) * t) * t) * t) *
               gammalib::deg2rad;

    // Compute Sun's mean anomaly in radians
    double M = (357.5291092 +
               (35999.0502909 +
               (-0.0001536 + 1.0/24490000.0 * t) * t) * t) *
               gammalib::deg2rad;

    // Compute Moon's mean anomaly in radians
    double Mprime = (134.9633964 +
                    (477198.8675055 +
                    (0.0087414 +
                    (1.0/69699.0 - 1.0/14712000.0 * t) * t) * t) * t) *
                    gammalib::deg2rad;

    // Compute Moon's argument of latitude in radians
    double F = (93.2720950 +
               (483202.0175233 +
               (-0.0036539 +
               (-1.0/35260000 + 1.0/863310000 * t) * t) * t) * t) *
               gammalib::deg2rad;

    // Compute longitude of the ascending node of the Moon's mean orbit on
    // the ecliptic, measured from the mean equinox of the date in radians
    double Omega = (125.04452 +
                   (-1934.136261 +
                   (0.0020708 + 1.0/4.5e5 * t) * t) * t) *
                   gammalib::deg2rad;

    // Compute actions of Venus and Jupiter in radians
    double A1 = (119.75 +    131.849 * t) * gammalib::deg2rad;
    double A2 = ( 53.09 + 479264.290 * t) * gammalib::deg2rad;
    double A3 = (313.45 + 481266.484 * t) * gammalib::deg2rad;

    // Compute sum quantities
    double Al =  3958.0 * std::sin(A1) +
                 1962.0 * std::sin(Lprime - F) +
                  318.0 * std::sin(A2);
    double Ab = -2235.0 * std::sin(Lprime) +
                  382.0 * std::sin(A3) +
                  175.0 * std::sin(A1 - F) +
                  175.0 * std::sin(A1 + F) +
                  127.0 * std::sin(Lprime - Mprime) -
                  115.0 * std::sin(Lprime + Mprime);

    // Compute eccentricity of Earth's orbit around the Sun
    double E  = 1.0 - (0.002516 + 7.4e-6 * t) * t;
    double E2 = E * E;

    // Sum the periodic terms
    double sum_l = Al;
    double sum_b = Ab;
    double sum_r = 0.0; // A_r = 0
    for (int i = 0; i < 60; ++i) {

        // Compute sine and cosine of longitude and sine of latitude
        double sinlng = sin_lng[i];
        //double coslng = cos_lng[i];
        double sinlat = sin_lat[i];
        if (std::abs(m_lng[i]) == 1) {
            sinlng *= E;
        //    coslng *= E;
        }
        else if (std::abs(m_lng[i]) == 2) {
            sinlng *= E2;
        //    coslng *= E2;
        }
        if (std::abs(m_lat[i]) == 1) {
            sinlat *= E;
        }
        else if (std::abs(m_lat[i]) == 2) {
            sinlat *= E2;
        }

        // Compute coefficient
        double s_l = d_lng[i] * D + m_lng[i] * M + mp_lng[i] * Mprime +
                     f_lng[i] * F;
        double s_b = d_lat[i] * D + m_lat[i] * M + mp_lat[i] * Mprime +
                     f_lat[i] * F;

        // Compute longitude
        sum_l += sinlng * std::sin(s_l);
        sum_b += sinlat * std::sin(s_b);
        //sum_r += coslng * std::sin(s_l);

    } // endfor: looped over sum terms

    // Compute Moon's ecliptic coordinates in degrees and distance in km
    double geolong  = Lprimed + 1.0e-6 * sum_l;
    double geolat   = 1.0e-6 * sum_b;
    //double distance = 385000.56 + 1.0e-3 * sum_r;

    // Compute nutation in arcsec, based on Chapter 22 of "Astronomical
    // Algorithms" by Jean Meeus (1998, 2nd ed.) which is based on the 1980
    // IAU Theory of Nutation and includes all terms larger than 0.0003 arcsec.
    double nut_long  = 0.0;
    double nut_obliq = 0.0;
    for (int i = 0; i < 63; ++i) {
        double arg = nutate_d_lng[i]  * D      + nutate_m_lng[i]  * M +
                     nutate_mp_lng[i] * Mprime + nutate_f_lng[i]  * F +
                     nutate_om_lng[i] * Omega;
        nut_long  += (nutate_sdelt[i] * t + nutate_sin_lng[i]) * std::sin(arg);
        nut_obliq += (nutate_cdelt[i] * t + nutate_cos_lng[i]) * std::cos(arg);
    }
    nut_long  *= 0.0001;
    nut_obliq *= 0.0001;

    // Add nutation in longitude
    geolong += nut_long / 3600.0;
    
    // nutate, jd, nlong, elong
    // geolong= geolong + nlong/3.6d3
    // cirrange,geolong
    double lambda = geolong * gammalib::deg2rad;
    double beta   = geolat  * gammalib::deg2rad;

    // Find mean obliquity using J. Laskar's formula. For the formula, see
    // http://www.neoprogrammics.com/obliquity_of_the_ecliptic/
    double y       = t / 100.0;    // Time in years
    double epsilon = (84381.448 +  // In degrees
                     (4680.93 +
                     (-1.55   +
                     (1999.25 +
                     (-51.38  +
                     (-249.67 +
                     (-39.05  +
                     (7.12    +
                     (27.87   +
                     (5.79 + 2.45 * y)*y)*y)*y)*y)*y)*y)*y)*y)*y)/3600.0;
    double eps     = (epsilon + nut_obliq/3600.0) * gammalib::deg2rad;

    // Convert (lambda,beta) in radians to (RA,Dec) in degrees
    double ra  = std::atan2(std::sin(lambda) * std::cos(eps) -
                            std::tan(beta)   * std::sin(eps), std::cos(lambda)) *
                 gammalib::rad2deg;
    double dec = std::asin(std::sin(beta)    * std::cos(eps) +
                           std::cos(beta)    * std::sin(eps) * std::sin(lambda)) *
                 gammalib::rad2deg;

    // Set position
    radec_deg(ra, dec);

    // Precess sky position to requested epoch
    precess(time.julian_epoch(), epoch);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return sky direction as 3D vector in celestial coordinates
 *
 * @return Sky direction as 3D vector in celestial coordinates.
 ***************************************************************************/
GVector GSkyDir::celvector(void) const
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb) {
        gal2equ();
    }

    // Compute 3D vector
    double  cosra  = std::cos(m_ra);
    double  sinra  = std::sin(m_ra);
    #if defined(G_SINCOS_CACHE)
    if (!m_has_radec_cache) {
        m_has_radec_cache = true;
        m_sin_dec         = std::sin(m_dec);
        m_cos_dec         = std::cos(m_dec);
    }
    GVector vector(m_cos_dec*cosra, m_cos_dec*sinra, m_sin_dec);
    #else
    double  cosdec = std::cos(m_dec);
    double  sindec = std::sin(m_dec);
    GVector vector(cosdec*cosra, cosdec*sinra, sindec);
    #endif

    // Return vector
    return vector;
}


/***********************************************************************//**
 * @brief Compute cosine of angular distance between sky directions
 *
 * @param[in] dir Sky direction to which cosine of distance is to be computed.
 * @return Cosine of angular distance.
 *
 * Computes the cosine of the angular distance between two sky directions in
 * radians.
 ***************************************************************************/
double GSkyDir::cos_dist(const GSkyDir& dir) const
{
    // Initialise cosine of distance
    double cosdis;

    // Compute dependent on coordinate system availability. This speeds
    // up things by avoiding unnecessary coordinate transformations.
    if (m_has_lb) {
        #if defined(G_SINCOS_CACHE)
        if (!m_has_lb_cache) {
            m_has_lb_cache = true;
            m_sin_b        = std::sin(m_b);
            m_cos_b        = std::cos(m_b);
        }
        #endif
        if (dir.m_has_lb) {
            #if defined(G_SINCOS_CACHE)
            if (!dir.m_has_lb_cache) {
                dir.m_has_lb_cache = true;
                dir.m_sin_b        = std::sin(dir.m_b);
                dir.m_cos_b        = std::cos(dir.m_b);
            }
            cosdis = m_sin_b * dir.m_sin_b +
                     m_cos_b * dir.m_cos_b *
                     std::cos(dir.m_l - m_l);
            #else
            cosdis = std::sin(m_b) * std::sin(dir.m_b) +
                     std::cos(m_b) * std::cos(dir.m_b) *
                     std::cos(dir.m_l - m_l);
            #endif
        }
        else {
            #if defined(G_SINCOS_CACHE)
            cosdis = m_sin_b * std::sin(dir.b()) +
                     m_cos_b * std::cos(dir.b()) *
                     std::cos(dir.l() - m_l);
            #else
            cosdis = std::sin(m_b) * std::sin(dir.b()) +
                     std::cos(m_b) * std::cos(dir.b()) *
                     std::cos(dir.l() - m_l);
            #endif
        }
    }
    else if (m_has_radec) {
        #if defined(G_SINCOS_CACHE)
        if (!m_has_radec_cache) {
            m_has_radec_cache = true;
            m_sin_dec         = std::sin(m_dec);
            m_cos_dec         = std::cos(m_dec);
        }
        #endif
        if (dir.m_has_radec) {
            #if defined(G_SINCOS_CACHE)
            if (!dir.m_has_radec_cache) {
                dir.m_has_radec_cache = true;
                dir.m_sin_dec         = std::sin(dir.m_dec);
                dir.m_cos_dec         = std::cos(dir.m_dec);
            }
            cosdis = m_sin_dec * dir.m_sin_dec +
                     m_cos_dec * dir.m_cos_dec *
                     std::cos(dir.m_ra - m_ra);
            #else
            cosdis = std::sin(m_dec) * std::sin(dir.m_dec) +
                     std::cos(m_dec) * std::cos(dir.m_dec) *
                     std::cos(dir.m_ra - m_ra);
            #endif
        }
        else {
            #if defined(G_SINCOS_CACHE)
            cosdis = m_sin_dec * sin(dir.dec()) +
                     m_cos_dec * cos(dir.dec()) *
                     std::cos(dir.ra() - m_ra);
            #else
            cosdis = sin(m_dec) * sin(dir.dec()) +
                     cos(m_dec) * cos(dir.dec()) *
                     std::cos(dir.ra() - m_ra);
            #endif
        }
    }
    else {
        cosdis = std::sin(dec()) * std::sin(dir.dec()) +
                 std::cos(dec()) * std::cos(dir.dec()) *
                 std::cos(dir.ra() - ra());
    }

    // Return cosine of distance
    return cosdis;
}


/***********************************************************************//**
 * @brief Compute position angle between sky directions in radians
 *
 * @param[in] dir Sky direction.
 * @return Position angle in radians.
 *
 * Computes the position angle using
 * \f[PA = \arctan \left( 
 *         \frac{\sin( \alpha_1 - \alpha_0 )}
 *              {\cos \delta_0 \tan \delta_1 -
 *               \sin \delta_0 \cos(\alpha_1 - \alpha_0)} \right)\f]
 * where
 * \f$(\alpha_0,\delta_0)\f$ are the coordinates of the reference point, and
 * \f$(\alpha_1,\delta_1)\f$ are the coordinates of sky direction for which
 * the position angle is to be computed.
 *
 * The position angle is counted counterclockwise from north.
 ***************************************************************************/
double GSkyDir::posang(const GSkyDir& dir) const
{
    // Initialise arguments of arctan
    double arg_1;
    double arg_2;

    // Compute dependent on coordinate system availability. This speeds
    // up things by avoiding unnecessary coordinate transformations.
    if (m_has_lb) {
        #if defined(G_SINCOS_CACHE)
        if (!m_has_lb_cache) {
            m_has_lb_cache = true;
            m_sin_b        = std::sin(m_b);
            m_cos_b        = std::cos(m_b);
        }
        #endif
        if (dir.m_has_lb) {
            arg_1 = std::sin(dir.m_l - m_l);
            #if defined(G_SINCOS_CACHE)
            arg_2 = m_cos_b * std::tan(dir.m_b) -
                    m_sin_b * std::cos(dir.m_l - m_l);
            #else
            arg_2 = std::cos(m_b) * std::tan(dir.m_b) -
                    std::sin(m_b) * std::cos(dir.m_l - m_l);
            #endif
        }
        else {
            arg_1 = std::sin(dir.l() - m_l);
            #if defined(G_SINCOS_CACHE)
            arg_2 = m_cos_b * std::tan(dir.b()) -
                    m_sin_b * std::cos(dir.l() - m_l);
            #else
            arg_2 = std::cos(m_b) * std::tan(dir.b()) -
                    std::sin(m_b) * std::cos(dir.l() - m_l);
            #endif
        }
    }
    else if (m_has_radec) {
        #if defined(G_SINCOS_CACHE)
        if (!m_has_radec_cache) {
            m_has_radec_cache = true;
            m_sin_dec         = std::sin(m_dec);
            m_cos_dec         = std::cos(m_dec);
        }
        #endif
        if (dir.m_has_radec) {
            arg_1 = std::sin(dir.m_ra - m_ra);
            #if defined(G_SINCOS_CACHE)
            arg_2 = m_cos_dec * std::tan(dir.m_dec) -
                    m_sin_dec * std::cos(dir.m_ra - m_ra);
            #else
            arg_2 = std::cos(m_dec) * std::tan(dir.m_dec) -
                    std::sin(m_dec) * std::cos(dir.m_ra - m_ra);
            #endif
        }
        else {
            arg_1 = std::sin(dir.ra() - m_ra);
            #if defined(G_SINCOS_CACHE)
            arg_2 = m_cos_dec * std::tan(dir.dec()) -
                    m_sin_dec * std::cos(dir.ra() - m_ra);
            #else
            arg_2 = std::cos(m_dec) * std::tan(dir.dec()) -
                    std::sin(m_dec) * std::cos(dir.ra() - m_ra);
            #endif
        }
    }
    else {
        arg_1 = std::sin(dir.ra() - ra());
        arg_2 = std::cos(dec())*std::tan(dir.dec()) -
                std::sin(dec())*std::cos(dir.ra() - ra());
    }

    // Compute position angle
    double pa = std::atan2(arg_1, arg_2);

    // Return position angle
    return pa;
}


/***********************************************************************//**
 * @brief Print sky direction information
 *
 * @param[in] chatter Chattiness.
 * @return String containing sky direction information.
 ***************************************************************************/
std::string GSkyDir::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Put coordinates in string
        if (m_has_lb) {
            result = "(l,b)=("+gammalib::str(m_l*gammalib::rad2deg) + "," +
                     gammalib::str(m_b*gammalib::rad2deg)+")";
        }
        else if (m_has_radec) {
            result = "(RA,Dec)=("+gammalib::str(m_ra*gammalib::rad2deg) + "," +
                     gammalib::str(m_dec*gammalib::rad2deg)+")";
        }
        else {
            result = "(RA,Dec)=(not initialised)";
        }

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
void GSkyDir::init_members(void)
{
    // Initialise members
    m_has_lb    = false;
    m_has_radec = false;
    m_l         = 0.0;
    m_b         = 0.0;
    m_ra        = 0.0;
    m_dec       = 0.0;

    // Initialise sincos cache
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = false;
    m_has_radec_cache = false;
    m_sin_b           = 0.0;
    m_cos_b           = 0.0;
    m_sin_dec         = 0.0;
    m_cos_dec         = 0.0;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Sky direction.
 ***************************************************************************/
void GSkyDir::copy_members(const GSkyDir& dir)
{
    // Copy members
    m_has_lb    = dir.m_has_lb;
    m_has_radec = dir.m_has_radec;
    m_l         = dir.m_l;
    m_b         = dir.m_b;
    m_ra        = dir.m_ra;
    m_dec       = dir.m_dec;

    // Copy sincos cache
    #if defined(G_SINCOS_CACHE)
    m_has_lb_cache    = dir.m_has_lb_cache;
    m_has_radec_cache = dir.m_has_radec_cache;
    m_sin_b           = dir.m_sin_b;
    m_cos_b           = dir.m_cos_b;
    m_sin_dec         = dir.m_sin_dec;
    m_cos_dec         = dir.m_cos_dec;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyDir::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert equatorial to galactic coordinates
 ***************************************************************************/
void GSkyDir::equ2gal(void) const
{
    // Get non-const pointers to data members. This allows to circumvent
    // the const correctness and allows treating GSkyDir access methods
    // as const
    double* l = (double*)&m_l;
    double* b = (double*)&m_b;

    // Convert from equatorial to galactic
    euler(0, m_ra, m_dec, l, b);

    // Signal that galactic coordinates are available
    m_has_lb = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert galactic to equatorial coordinates
 ***************************************************************************/
void GSkyDir::gal2equ(void) const
{
    // Get non-const pointers to data members. This allows to circumvent
    // the const correctness and allows treating GSkyDir access methods
    // as const
    double* ra  = (double*)&m_ra;
    double* dec = (double*)&m_dec;

    // Convert from galactic to equatorial
    euler(1, m_l, m_b, ra, dec);

    // Signal that equatorial coordinates are available
    m_has_radec = true;

    // Return
    return;
}


/***********************************************************************//**
 * @brief General coordinate transformation routine for J2000
 *
 * @param[in] type Conversion type (0=equ2gal, 1=gal2equ)
 * @param[in] xin Input longitude (RA or GLON) in radians.
 * @param[in] yin Input latitude (Dec or GLAT) in radians.
 * @param[out] xout Output longitude in radians.
 * @param[out] yout Output latitude in radians.
 ***************************************************************************/
void GSkyDir::euler(const int& type, const double& xin, const double &yin, 
                    double* xout, double *yout) const
{
    // Set transformation constants
    const double psi[]    = {0.57477043300,  4.9368292465};
    const double stheta[] = {0.88998808748, -0.88998808748};
    const double ctheta[] = {0.45598377618,  0.45598377618};
    const double phi[]    = {4.9368292465,   0.57477043300};

    // Perform transformation
    double a    = xin - phi[type];
    double b    = yin;
    double sb   = std::sin(b);
    double cb   = std::cos(b);
    double cbsa = cb * std::sin(a);

    //
    a = std::atan2(ctheta[type] * cbsa + stheta[type] * sb, cb * std::cos(a));
    b = -stheta[type] * cbsa + ctheta[type] * sb;
    if (b > 1.0)
        b = 1.0;

    //
    *yout = std::asin(b);
    *xout = gammalib::modulo((a+psi[type] + gammalib::fourpi), gammalib::twopi);

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
 * @param[in] a First sky direction.
 * @param[in] b Second sky direction.
 *
 * Compare two sky directions. If the coordinate is at the pole, the Right
 * Ascension or Longitude value is irrelevant.
 *
 * Comparisons are done dependent on the available coordinate system. This
 * speeds up things and avoids unnecessary coordinate transformations.
 ***************************************************************************/
bool operator==(const GSkyDir &a, const GSkyDir &b)
{
    // Initialise result
    bool equal = false;

    // Compute dependent on coordinate system availability. This speeds
    // up things by avoiding unnecessary coordinate transformations.

    // Check if both have equatorial coordinates
    if (a.m_has_radec && b.m_has_radec) {
        if (std::abs(a.m_dec) == 90.0) {
            equal = (a.m_dec == b.m_dec);
        }
        else {
            equal = (a.m_dec == b.m_dec && a.m_ra == b.m_ra);
        }
    }
    // ... check if both have Galactic coordinates
    else if (a.m_has_lb && b.m_has_lb) {
        if (std::abs(a.m_b) == 90.0) {
            equal = (a.m_b == b.m_b);
        }
        else {
            equal = (a.m_b == b.m_b && a.m_l == b.m_l);
        }
    }
    // ... otherwise the coordinate systems are different
    else {
        if (a.m_has_lb) {
            if (std::abs(a.m_b) == 90.0) {
                equal = (a.m_b == b.b());
            }
            else {
                equal = (a.m_b == b.b() && a.m_l == b.l());
            }
        } 
        else if (a.m_has_radec) {
            if (std::abs(a.m_dec) == 90.0) {
                equal = (a.m_dec == b.dec());
            }
            else {
                equal = (a.m_dec == b.dec() && a.m_ra == b.ra());
            }
        }
        else {
            if (std::abs(b.dec()) == 90.0) {
                equal = (b.dec() == a.dec());
            }
            else {
                equal = (b.dec() == a.dec() && b.ra() == a.ra());
            }
        }
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
bool operator!=(const GSkyDir &a, const GSkyDir &b)
{
    // Return non equality
    return (!(a==b));
}
