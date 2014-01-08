/***************************************************************************
 *          GSkyProjection.cpp - Abstract sky projection base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GSkyProjection.cpp
 * @brief Abstract sky projection base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GSkyProjection.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COORDSYS_SET                "GSkyProjection::coordsys(std::string)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Local prototypes ___________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GSkyProjection::GSkyProjection(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] proj Sky projection.
 ***************************************************************************/
GSkyProjection::GSkyProjection(const GSkyProjection& proj)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(proj);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyProjection::~GSkyProjection(void)
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
 * @param[in] proj Sky projection.
 * @return Sky projection.
 ***************************************************************************/
GSkyProjection& GSkyProjection::operator= (const GSkyProjection& proj)
{
    // Execute only if object is not identical
    if (this != &proj) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(proj);

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
 * @brief Returns coordinate system.
 *
 * @return Coordinate system string.
 *
 * Returns one of 
 * 'EQU' (equatorial),
 * 'GAL' (galactic),
 ***************************************************************************/
std::string GSkyProjection::coordsys(void) const
{
    // Set coordinate system
    std::string s_coordsys;
    switch (m_coordsys) {
    case 0:
        s_coordsys = "EQU";
        break;
    case 1:
        s_coordsys = "GAL";
        break;
    default:
        s_coordsys = "UNKNOWN";
        break;
    }

    // Return coordinate system
    return s_coordsys;
}


/***********************************************************************//**
 * @brief Set coordinate system
 *
 * @param[in] coordsys Coordinate system
 *
 * @exception GException::wcs_bad_coords
 *            Invalid coordsys parameter.
 *
 * Set coordinate system from std::string. The method recognizes the
 * following codes:
 * 'EQU', 'CEL', 'C': celestial,
 * 'GAL', 'G': galactic,
 ***************************************************************************/
void GSkyProjection::coordsys(const std::string& coordsys)
{
    // Convert argument to upper case
    std::string ucoordsys = gammalib::toupper(coordsys);

    // Set coordinate system
    if (ucoordsys == "EQU" || ucoordsys == "CEL" || ucoordsys == "C") {
        m_coordsys = 0;
    }
    else if (ucoordsys == "GAL" || ucoordsys == "G") {
        m_coordsys = 1;
    }
    else {
        throw GException::wcs_bad_coords(G_COORDSYS_SET, coordsys);
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GSkyProjection::init_members(void)
{
    // Initialise members
    m_coordsys = 0; // 0 means EQU

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] proj Sky projection.
 ***************************************************************************/
void GSkyProjection::copy_members(const GSkyProjection& proj)
{
    // Copy attributes
    m_coordsys = proj.m_coordsys;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyProjection::free_members(void)
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
 * @param[in] a First sky projection.
 * @param[in] b Second sky projection.
 * @return True if @p a and @p b are identical.
 ***************************************************************************/
bool operator==(const GSkyProjection &a, const GSkyProjection &b)
{
    // Return result
    return a.compare(b);
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] a First sky projection.
 * @param[in] b Second sky projection.
 * @return True if @p a and @p b are not identical.
 ***************************************************************************/
bool operator!=(const GSkyProjection &a, const GSkyProjection &b)
{
    // Return result
    return !(a == b);
}
