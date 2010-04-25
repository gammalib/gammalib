/***************************************************************************
 *          GSkyDir.cpp  -  Class that implements a sky direction          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GSkyDir.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototype __________________________________________________________ */
void   euler(const int& type, const double& xin, const double &yin, 
             double* xout, double *yout);

/*==========================================================================
 =                                                                         =
 =                      GSkyDir constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GSkyDir::GSkyDir()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Sky direction from which class should be instantiated.
 ***************************************************************************/
GSkyDir::GSkyDir(const GSkyDir& dir)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyDir::~GSkyDir()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GSkyDir operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 ***************************************************************************/
GSkyDir& GSkyDir::operator= (const GSkyDir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(dir);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                          GSkyDir public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set equatorial sky direction (radians)
 *
 * @param[in] ra Right Ascension in radians.
 * @param[in] dec Declination in radians.
 ***************************************************************************/
void GSkyDir::radec(const double& ra, const double& dec)
{
    // Set attributes
    m_has_lb    = 0;
    m_has_radec = 1;

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
    m_has_lb    = 0;
    m_has_radec = 1;

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
    m_has_lb    = 1;
    m_has_radec = 0;

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
    m_has_lb    = 1;
    m_has_radec = 0;

    // Set direction
    m_l = l * deg2rad;
    m_b = b * deg2rad;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns galactic longitude in radians
 ***************************************************************************/
double GSkyDir::l(void)
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
double GSkyDir::l_deg(void)
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
double GSkyDir::b(void)
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
double GSkyDir::b_deg(void)
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
double GSkyDir::ra(void)
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
double GSkyDir::ra_deg(void)
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
double GSkyDir::dec(void)
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
double GSkyDir::dec_deg(void)
{
    // If we have no equatorial coordinates then get them now
    if (!m_has_radec && m_has_lb)
        gal2equ();
    
    // Return Declination
    return m_dec * rad2deg;
}


/*==========================================================================
 =                                                                         =
 =                          GSkyDir private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyDir::init_members(void)
{
    // Initialise members
    m_has_lb    = 0;
    m_has_radec = 0;
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
 * @param[in] dir Sky direction from which members should be copied
 ***************************************************************************/
void GSkyDir::copy_members(const GSkyDir& dir)
{
    // Copy attributes
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
    // Free memory

    // Signal free pointers

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert equatorial to galactic coordinates
 ***************************************************************************/
void GSkyDir::equ2gal(void)
{
    // Convert from equatorial to galactic
    euler(0, m_ra, m_dec, &m_l, &m_b);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert galactic to equatorial coordinates
 ***************************************************************************/
void GSkyDir::gal2equ(void)
{
    // Convert from galactic to equatorial
    euler(1, m_l, m_b, &m_ra, &m_dec);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              GSkyDir friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] column Sky direction to put in output stream
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSkyDir& dir)
{
    // Create coordinate system dependent output
    if (dir.m_has_lb)
        os << "(l,b)=(" << dir.m_l*rad2deg << "," << dir.m_b*rad2deg << ")";
    else if (dir.m_has_radec)
        os << "(RA,Dec)=(" << dir.m_ra*rad2deg << "," << dir.m_dec*rad2deg << ")";

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GSkyDir                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief General coordinate transformation routine for J2000
 *
 * @param[in] type Conversion type (0=equ2gal, 1=gal2equ)
 * @param[in] xin Input longitude (RA or GLON) in radians.
 * @param[in] yin Input latitude (Dec or GLAT) in radians.
 * @param[out] xout Output longitude in radians.
 * @param[out] yout Output latitude in radians.
 ***************************************************************************/
void euler(const int& type, const double& xin, const double &yin, 
           double* xout, double *yout)
{
    // Set transformation constants
    const double psi[]    = {0.57477043300,  4.9368292465};
    const double stheta[] = {0.88998808748, -0.88998808748};
    const double ctheta[] = {0.45598377618,  0.45598377618};
    const double phi[]    = {4.9368292465,   0.57477043300};
    
    // Perform transformation
    double a    = xin - phi[type];
    double b    = yin;
    double sb   = sin(b);
    double cb   = cos(b);
    double cbsa = cb * sin(a);
    
    //
    a = atan2(ctheta[type] * cbsa + stheta[type] * sb, cb * cos(a));
    b = -stheta[type] * cbsa + ctheta[type] * sb;
    if (b > 1.0)
        b = 1.0;

    //
    *yout = asin(b);
    *xout = modulo((a+psi[type] + fourpi), twopi);

    // Return
    return;
}
