/***************************************************************************
 *    GModelSpectralTablePar.cpp - Spectral table model parameter class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Juergen Knoedlseder                         *
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
 * @file GModelSpectralTablePar.cpp
 * @brief Spectral table model parameter class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <algorithm>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelPar.hpp"
#include "GModelSpectralTablePar.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */

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
GModelSpectralTablePar::GModelSpectralTablePar(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Model parameter constructor
 *
 * @param[in] par Model parameter.
 * @param[in] values Parameter values.
 *
 * Constructs a model table parameter combining a model parameter with a
 * vector of parameter values. The values in the @p values vector may be
 * unsorted, the constructor will put the values into an acsending order.
 ***************************************************************************/
GModelSpectralTablePar::GModelSpectralTablePar(const GModelPar&           par,
                                               const std::vector<double>& values)
{
    // Initialise members
    init_members();

    // Sort parameter values
    std::vector<double> sorted_values = values;
    std::sort(sorted_values.begin(), sorted_values.end());

    // Set members
    m_par    = par;
    m_values = GNodeArray(sorted_values);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] par Table model parameter.
 ***************************************************************************/
GModelSpectralTablePar::GModelSpectralTablePar(const GModelSpectralTablePar& par)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpectralTablePar::~GModelSpectralTablePar(void)
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
 * @param[in] par Table model parameter.
 * @return Table model parameter.
 ***************************************************************************/
GModelSpectralTablePar& GModelSpectralTablePar::operator=(const GModelSpectralTablePar& par)
{
    // Execute only if object is not identical
    if (this != &par) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(par);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear table model parameter
***************************************************************************/
void GModelSpectralTablePar::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone table model parameter
***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePar::clone(void) const
{
    // Clone table model parameter
    return new GModelSpectralTablePar(*this);
}


/***********************************************************************//**
 * @brief Print table model parameter
 *
 * @param[in] chatter Chattiness.
 * @return String containing table model parameter information.
 ***************************************************************************/
std::string GModelSpectralTablePar::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralTablePar ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name"));
        result.append(m_par.name());
        result.append("\n"+gammalib::parformat("Number of values"));
        result.append(gammalib::str(m_values.size()));

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
void GModelSpectralTablePar::init_members(void)
{
    // Initialize members
    m_par.clear();
    m_values.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] par Table model parameter.
 ***************************************************************************/
void GModelSpectralTablePar::copy_members(const GModelSpectralTablePar& par)
{
    // Copy members
    m_par    = par.m_par;
    m_values = par.m_values;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralTablePar::free_members(void)
{
    // Return
    return;
}
