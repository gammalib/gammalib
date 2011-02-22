/***************************************************************************
 *          GSkyDir.cpp  -  Class that implements a sky direction          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkyDir.cpp
 * @brief Sky direction class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GSkyDir.hpp"
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
 * @brief Set equatorial sky direction (radians)
 *
 * @param[in] ra Right Ascension in radians.
 * @param[in] dec Declination in radians.
 ***************************************************************************/
void GSkyDir::radec(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = false;
    m_has_radec = true;

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
 ***************************************************************************/
void GSkyDir::radec_deg(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = false;
    m_has_radec = true;

    // Set direction
    m_ra  = ra  * deg2rad;
    m_dec = dec * deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set galactic sky direction (radians)
 *
 * @param[in] l Galactic longitude in radians.
 * @param[in] b Galactic latitude in radians.
 ***************************************************************************/
void GSkyDir::lb(const double& l, const double& b)
{
    // Set attributes
    m_has_lb    = true;
    m_has_radec = false;

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
 ***************************************************************************/
void GSkyDir::lb_deg(const double& l, const double& b)
{
    // Set attributes
    m_has_lb    = true;
    m_has_radec = false;

    // Set direction
    m_l = l * deg2rad;
    m_b = b * deg2rad;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sky direction from 3D vector in celestial coordinates
 *
 * @param[in] vector 3D vector.
 ***************************************************************************/
void GSkyDir::celvector(const GVector& vector)
{
    // Set attributes
    m_has_lb    = false;
    m_has_radec = true;

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
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Allocate Euler and rotation matrices
    GMatrix ry;
    GMatrix rz;
    GMatrix rot;

    // Set up rotation matrix to rotate from native coordinates to
    // celestial coordinates
    ry.eulery(m_dec*rad2deg - 90.0);
    rz.eulerz(-m_ra*rad2deg);
    rot = transpose(ry * rz);

    // Set up native coordinate vector
    double phi_rad   = phi   * deg2rad;
    double theta_rad = theta * deg2rad;
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
 * @brief Returns galactic longitude in radians
 ***************************************************************************/
double GSkyDir::l(void) const
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic longitude
    return m_l;
}


/***********************************************************************//**
 * @brief Returns galactic longitude in degrees
 ***************************************************************************/
double GSkyDir::l_deg(void) const
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic longitude
    return m_l * rad2deg;
}


/***********************************************************************//**
 * @brief Returns galactic latitude in radians
 ***************************************************************************/
double GSkyDir::b(void) const
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic latitude
    return m_b;
}


/***********************************************************************//**
 * @brief Returns galactic latitude in degrees
 ***************************************************************************/
double GSkyDir::b_deg(void) const
{
    // If we have no galactic coordinates then get them now
    if (!m_has_lb && m_has_radec)
        equ2gal();

    // Return galactic latitude
    return m_b * rad2deg;
}


/***********************************************************************//**
 * @brief Returns Right Ascension in radians
 ***************************************************************************/
double GSkyDir::ra(void) const
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Right Ascension
    return m_ra;
}


/***********************************************************************//**
 * @brief Returns Right Ascension in degrees
 ***************************************************************************/
double GSkyDir::ra_deg(void) const
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Right Ascension
    return m_ra * rad2deg;
}


/***********************************************************************//**
 * @brief Returns Declination in radians
 ***************************************************************************/
double GSkyDir::dec(void) const
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Declination
    return m_dec;
}


/***********************************************************************//**
 * @brief Returns Declination in degrees
 ***************************************************************************/
double GSkyDir::dec_deg(void) const
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Return Declination
    return m_dec * rad2deg;
}


/***********************************************************************//**
 * @brief Returns sky direction as 3D vector in celestial coordinates
 ***************************************************************************/
GVector GSkyDir::celvector(void) const
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();

    // Compute 3D vector
    double  cosra  = std::cos(m_ra);
    double  sinra  = std::sin(m_ra);
    double  cosdec = std::cos(m_dec);
    double  sindec = std::sin(m_dec);
    GVector vector(cosdec*cosra, cosdec*sinra, sindec);

    // Return vector
    return vector;
}


/***********************************************************************//**
 * @brief Compute angular distance between sky directions in radians
 *
 * @param[in] dir Sky direction to which distance is to be computed.
 ***************************************************************************/
double GSkyDir::dist(const GSkyDir& dir) const
{
    // Initialise cosine of distance
    double cosdis;

    // Compute dependent on coordinate system availability. This speeds
    // up things by avoiding unnecessary coordinate transformations.
    if (m_has_lb && dir.m_has_lb) {
        cosdis = std::sin(m_b)*std::sin(dir.m_b) +
                 std::cos(m_b)*std::cos(dir.m_b) * std::cos(dir.m_l - m_l);
    }
    else if (m_has_radec && dir.m_has_radec) {
        cosdis = std::sin(m_dec)*std::sin(dir.m_dec) +
                 std::cos(m_dec)*std::cos(dir.m_dec) * std::cos(dir.m_ra - m_ra);
    }
    else if (m_has_lb) {
        cosdis = std::sin(m_b)*std::sin(dir.b()) +
                 std::cos(m_b)*std::cos(dir.b()) * std::cos(dir.l() - m_l);
    }
    else if (m_has_radec) {
        cosdis = sin(m_dec)*sin(dir.dec()) +
                 cos(m_dec)*cos(dir.dec()) * std::cos(dir.ra() - m_ra);
    }
    else {
        cosdis = std::sin(dec())*std::sin(dir.dec()) +
                 std::cos(dec())*std::cos(dir.dec()) * std::cos(dir.ra() - ra());
    }

    // Compute distance (use argument save GTools function)
    double dist = arccos(cosdis);

    // Return distance
    return dist;
}


/***********************************************************************//**
 * @brief Compute angular distance between sky directions in degrees
 *
 * @param[in] dir Sky direction to which distance is to be computed.
 ***************************************************************************/
double GSkyDir::dist_deg(const GSkyDir& dir) const
{
    // Return distance in degrees
    return (dist(dir) * rad2deg);
}


/***********************************************************************//**
 * @brief Compute position angle between sky directions in radians
 *
 * @param[in] dir Sky direction.
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
    if (m_has_lb && dir.m_has_lb) {
        arg_1 = std::sin(dir.m_l - m_l);
        arg_2 = std::cos(m_b)*std::tan(dir.m_b) - std::sin(m_b)*std::cos(dir.m_l - m_l);
    }
    else if (m_has_radec && dir.m_has_radec) {
        arg_1 = std::sin(dir.m_ra - m_ra);
        arg_2 = std::cos(m_dec)*std::tan(dir.m_dec) - std::sin(m_dec)*std::cos(dir.m_ra - m_ra);
    }
    else if (m_has_lb) {
        arg_1 = std::sin(dir.l() - m_l);
        arg_2 = std::cos(m_b)*std::tan(dir.b()) - std::sin(m_b)*std::cos(dir.l() - m_l);
    }
    else if (m_has_radec) {
        arg_1 = std::sin(dir.ra() - m_ra);
        arg_2 = std::cos(m_dec)*std::tan(dir.dec()) - std::sin(m_dec)*std::cos(dir.ra() - m_ra);
    }
    else {
        arg_1 = std::sin(dir.ra() - ra());
        arg_2 = std::cos(dec())*std::tan(dir.dec()) - std::sin(dec())*std::cos(dir.ra() - ra());
    }

    // Compute position angle
    double pa = std::atan2(arg_1, arg_2);

    // Return position angle
    return pa;
}


/***********************************************************************//**
 * @brief Compute position angle between sky directions in degrees
 *
 * @param[in] dir Sky direction.
 *
 * See GSkyDir::posang for more information about the computation of the
 * position angle.
 ***************************************************************************/
double GSkyDir::posang_deg(const GSkyDir& dir) const
{
    // Return position angle in degrees
    return (posang(dir) * rad2deg);
}


/***********************************************************************//**
 * @brief Print vector information
 ***************************************************************************/
std::string GSkyDir::print(void) const
{
    // Initialise result string
    std::string result;

    // Put coordinates in string
    if (m_has_lb)
        result = "(l,b)=("+str(m_l*rad2deg)+","+str(m_b*rad2deg)+")";
    else if (m_has_radec)
        result = "(RA,Dec)=("+str(m_ra*rad2deg)+","+str(m_dec*rad2deg)+")";
    else
        result = "(RA,Dec)=(not initialised)";

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
    *xout = modulo((a+psi[type] + fourpi), twopi);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] dir Sky direction.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSkyDir& dir)
{
     // Write sky direction in output stream
    os << dir.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] dir Sky direction.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GSkyDir& dir)
{
    // Write sky direction into logger
    log << dir.print();

    // Return logger
    return log;
}


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
    if (a.m_has_lb && b.m_has_lb) {
        if (std::abs(a.m_b) == 90.0)
            equal = (a.m_b == b.m_b);
        else
            equal = (a.m_b == b.m_b && a.m_l == b.m_l);
    }
    else if (a.m_has_radec && b.m_has_radec) {
        if (std::abs(a.m_dec) == 90.0)
            equal = (a.m_dec == b.m_dec);
        else
            equal = (a.m_dec == b.m_dec && a.m_ra == b.m_ra);
    }
    else if (a.m_has_lb) {
        if (std::abs(a.m_b) == 90.0)
            equal = (a.m_b == b.b());
        else
            equal = (a.m_b == b.b() && a.m_l == b.l());
    }
    else if (a.m_has_radec) {
        if (std::abs(a.m_dec) == 90.0)
            equal = (a.m_dec == b.dec());
        else
            equal = (a.m_dec == b.dec() && a.m_ra == b.ra());
    }
    else {
        if (std::abs(b.dec()) == 90.0)
            equal = (b.dec() == a.dec());
        else
            equal = (b.dec() == a.dec() && b.ra() == a.ra());
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

