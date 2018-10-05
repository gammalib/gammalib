/***************************************************************************
 *         GCTAModelRadial.cpp  -  Abstract radial model base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadial.cpp
 * @brief Abstract radial background model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GCTAModelRadial.hpp"
#include "GCTAInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                           "GModelSpatial::operator[](int&)"
#define G_ACCESS2                   "GModelSpatial::operator[](std::string&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelRadial::GCTAModelRadial(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Radial background model.
 ***************************************************************************/
GCTAModelRadial::GCTAModelRadial(const GCTAModelRadial& model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelRadial::~GCTAModelRadial(void)
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
 * @param[in] model Radial background model.
 ***************************************************************************/
GCTAModelRadial& GCTAModelRadial::operator=(const GCTAModelRadial& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] dir Event direction.
 * @param[in] energy Event energy (not used).
 * @param[in] time Event time (not used).
 * @param[in] gradients Compute gradients?
 * @return Function value
 *
 * Evaluate radial model for a given event direction. The energy and time of
 * the event are not used.
 ***************************************************************************/
double GCTAModelRadial::eval(const GCTAInstDir& dir,
                             const GEnergy&     energy,
                             const GTime&       time,
                             const bool&        gradients) const
{
    // Compute offset angle in degrees
    double offset = dir.theta() * gammalib::rad2deg;

    // Evaluate function
    double value = eval(offset, gradients);

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] energy Event energy (not used).
 * @param[in] time Event time (not used).
 * @param[in,out] ran Random number generator.
 * @return Instrument direction
 *
 * Return random instrument direction. The energy and time of the event are
 * not used.
 ***************************************************************************/
GCTAInstDir GCTAModelRadial::mc(const GEnergy& energy,
                                const GTime&   time,
                                GRan&          ran) const
{
    // Get random instrument direction
    GCTAInstDir dir = mc(ran);

    // Return instrument direction
    return dir;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelRadial::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
void GCTAModelRadial::copy_members(const GCTAModelRadial& model)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRadial::free_members(void)
{
    // Return
    return;
}
