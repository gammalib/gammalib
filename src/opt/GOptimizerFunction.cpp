/***************************************************************************
 *     GOptimizerFunction.cpp - Optimizer function abstract base class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GOptimizerFunction.cpp
 * @brief Optimizer function abstract base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GOptimizerFunction.hpp"

/* __ Method name definitions ____________________________________________ */

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
GOptimizerFunction::GOptimizerFunction(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] fct Optimizer function.
 ***************************************************************************/
GOptimizerFunction::GOptimizerFunction(const GOptimizerFunction& fct)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(fct);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GOptimizerFunction::~GOptimizerFunction()
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
 * @param[in] fct Optimizer function.
 * @return Optimizer function.
 ***************************************************************************/
GOptimizerFunction& GOptimizerFunction::operator=(const GOptimizerFunction& fct)
{
    // Execute only if object is not identical
    if (this != &fct) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(fct);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GOptimizerFunction::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fct Optimizer function.
 ***************************************************************************/
void GOptimizerFunction::copy_members(const GOptimizerFunction& fct)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerFunction::free_members(void)
{
    // Return
    return;
}
