/***************************************************************************
 *           GWcs.cpp  -  World Coordinate System virtual base class       *
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
 * @file GWcs.cpp
 * @brief World Coordinate System virtual base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GWcs.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COORDSYS_SET                          "GWcs::coordsys(std::string)"

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
GWcs::GWcs(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param wcs World Coordinate System.
 ***************************************************************************/
GWcs::GWcs(const GWcs& wcs)
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
GWcs::~GWcs(void)
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
GWcs& GWcs::operator= (const GWcs& wcs)
{
    // Execute only if object is not identical
    if (this != &wcs) {

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
 * @brief Returns coordinate system.
 *
 * Returns one of 
 * 'EQU' (equatorial),
 * 'GAL' (galactic),
 * 'ECL' (ecliptic),
 * 'HEL' (helioecliptic), and
 * 'SGL' (supergalactic).
 ***************************************************************************/
std::string GWcs::coordsys(void) const
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
    case 2:
        s_coordsys = "ECL";
        break;
    case 3:
        s_coordsys = "HEL";
        break;
    case 4:
        s_coordsys = "SGL";
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
 * 'ECL', 'E': ecliptic,
 * 'HEL', 'H': helioecliptic, and
 * 'SGL', 'S': helioecliptic.
 ***************************************************************************/
void GWcs::coordsys(const std::string& coordsys)
{
    // Convert argument to upper case
    std::string ucoordsys = toupper(coordsys);

    // Set coordinate system
    if (ucoordsys == "EQU" || ucoordsys == "CEL" || ucoordsys == "C")
        m_coordsys = 0;
    else if (ucoordsys == "GAL" || ucoordsys == "G")
        m_coordsys = 1;
    else if (ucoordsys == "ECL" || ucoordsys == "E")
        m_coordsys = 2;
    else if (ucoordsys == "HEL" || ucoordsys == "H")
        m_coordsys = 3;
    else if (ucoordsys == "SGL" || ucoordsys == "S")
        m_coordsys = 4;
    else
        throw GException::wcs_bad_coords(G_COORDSYS_SET, coordsys);

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
void GWcs::init_members(void)
{
    // Initialise members
    m_coordsys = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
void GWcs::copy_members(const GWcs& wcs)
{
    // Copy attributes
    m_coordsys = wcs.m_coordsys;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcs::free_members(void)
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
 * @param[in] a First World Coordinate System.
 * @param[in] b Second World Coordinate System.
 ***************************************************************************/
bool operator== (const GWcs &a, const GWcs &b)
{
    // Return result
    return a.compare(b);
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] a First World Coordinate System.
 * @param[in] b Second World Coordinate System.
 ***************************************************************************/
bool operator!= (const GWcs &a, const GWcs &b)
{
    // Return result
    return !(a == b);
}
