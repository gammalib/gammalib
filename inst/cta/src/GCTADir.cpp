/***************************************************************************
 *                   GCTADir.cpp  -  CTA direction class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTADir.cpp
 * @brief CTA camera direction class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GCTADir.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Prototypes _________________________________________________________ */

/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTADir::GCTADir(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Camera direction.
 ***************************************************************************/
GCTADir::GCTADir(const GCTADir& dir)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky direction constructor
 *
 * @param[in] dir Sky direction.
 * @param[in] pnt CTA pointing.
 ***************************************************************************/
GCTADir::GCTADir(const GSkyDir& dir, const GCTAPointing& pnt)
{
    // Initialise class members
    init_members();

    // Set members
    this->dir(dir, pnt);

    // Return
    return;
}


/***********************************************************************//**
 * @brief CTA instrument direction constructor
 *
 * @param[in] dir CTA instrument direction.
 * @param[in] pnt CTA pointing.
 ***************************************************************************/
GCTADir::GCTADir(const GCTAInstDir& dir, const GCTAPointing& pnt)
{
    // Initialise class members
    init_members();

    // Set members
    this->dir(dir, pnt);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTADir::~GCTADir(void)
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
 * @param[in] dir Camera direction.
 ***************************************************************************/
GCTADir& GCTADir::operator= (const GCTADir& dir)
{
    // Execute only if object is not identical
    if (this != &dir) {

        // Free members
        free_members();

        // Initialise private members
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
 * @brief Clear instance
 ***************************************************************************/
void GCTADir::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCTADir* GCTADir::clone(void) const
{
    return new GCTADir(*this);
}


/***********************************************************************//**
 * @brief Convert sky direction into camera direction
 *
 * @param[in] dir Sky direction.
 * @param[in] pnt CTA pointing.
 ***************************************************************************/
void GCTADir::dir(const GSkyDir& dir, const GCTAPointing& pnt)
{
    // Compute radian offset and polar angles
    m_theta = pnt.dir().dist(dir);
    m_phi   = pnt.dir().posang(dir);

    // Flag cache as non-valid
    m_has_cache = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert CTA instrument direction into camera direction
 *
 * @param[in] dir CTA instrument direction.
 * @param[in] pnt CTA pointing.
 ***************************************************************************/
void GCTADir::dir(const GCTAInstDir& dir, const GCTAPointing& pnt)
{
    // Use sky direction method
    this->dir(dir.skydir(), pnt);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return radial offset angle in degrees
 ***************************************************************************/
double GCTADir::theta_deg(void) const
{
    // Return value
    return (m_theta*rad2deg);
}


/***********************************************************************//**
 * @brief Return polar angle in degrees
 ***************************************************************************/
double GCTADir::phi_deg(void) const
{
    // Return value
    return (m_phi*rad2deg);
}


/***********************************************************************//**
 * @brief Return cosine of radial offset angle
 ***************************************************************************/
const double& GCTADir::costheta(void) const
{
    // Update cache
    update();

    // Return value
    return m_cos_theta;
}


/***********************************************************************//**
 * @brief Return sine of radial offset angle
 ***************************************************************************/
const double& GCTADir::sintheta(void) const
{
    // Update cache
    update();

    // Return value
    return m_sin_theta;
}


/***********************************************************************//**
 * @brief Print camera direction information
 ***************************************************************************/
std::string GCTADir::print(void) const
{
    // Initialise result string
    std::string result;

    // Append camera direction
    result.append("(theta,phi)=(");
    result.append(str(theta_deg())+","+str(phi_deg()));
    result.append(")");

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
void GCTADir::init_members(void)
{
    // Initialise members
    m_theta = 0.0;
    m_phi   = 0.0;

    // Initialise cache
    m_has_cache = false;
    m_cos_theta = 0.0;
    m_sin_theta = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dir Camera direction.
 ***************************************************************************/
void GCTADir::copy_members(const GCTADir& dir)
{
    // Copy attributes
    m_theta = dir.m_theta;
    m_phi   = dir.m_phi;

    // Copy cache
    m_has_cache = dir.m_has_cache;
    m_cos_theta = dir.m_cos_theta;
    m_sin_theta = dir.m_sin_theta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTADir::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update precomputation cache
 ***************************************************************************/
void GCTADir::update(void) const
{
    // If cache is not valid them perform precomputations
    if (!m_has_cache) {

        // Perform precomputations
        m_cos_theta = std::cos(m_theta);
        m_sin_theta = std::sin(m_theta);

        // Signal cache validity
        m_has_cache = true;

    } // endif: Cache update was required

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
 * @param[in] dir Camera direction.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCTADir& dir)
{
    // Write camera direction in output stream
    os << dir.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] dir Camera direction.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GCTADir& dir)
{
    // Write camera direction into logger
    log << dir.print();

    // Return logger
    return log;
}
