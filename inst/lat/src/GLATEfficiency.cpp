/***************************************************************************
 *      GLATEfficiency.cpp - Fermi-LAT IRF efficiency factor functor       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GLATEfficiency.cpp
 * @brief Fermi-LAT IRF efficiency factor functor class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATEfficiency.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATEfficiency::GLATEfficiency(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Assignment constructor
 *
 * @param[in] pars Efficiency factor parameters.
 ***************************************************************************/
GLATEfficiency::GLATEfficiency(const std::vector<double>& pars)
{
    // Initialise class members
    init_members();

    // Assign members
    m_a0     = pars.at(0);
    m_b0     = pars.at(1);
    m_a1     = pars.at(2);
    m_logEb1 = pars.at(3);
    m_a2     = pars.at(4);
    m_logEb2 = pars.at(5);

    // Precompute additional offsets
    m_b1 = (m_a0 - m_a1)*m_logEb1 + m_b0;
    m_b2 = (m_a1 - m_a2)*m_logEb2 + m_b1;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] eff Efficiency factor functor.
 ***************************************************************************/
GLATEfficiency::GLATEfficiency(const GLATEfficiency& eff)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(eff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATEfficiency::~GLATEfficiency(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] eff Efficiency factor functor.
 * @return Efficiency factor functor.
 ***************************************************************************/
GLATEfficiency& GLATEfficiency::operator=(const GLATEfficiency& eff)
{
    // Execute only if object is not identical
    if (this != &eff) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(eff);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Functor operator
 *
 * @param[in] logE log10 of energy (MeV).
 ***************************************************************************/
double GLATEfficiency::operator() (const double& logE) const
{
    // Initialise the result
    double factor = 0;
    
    // Set efficiency factor depending on the energy domain
    if (logE < m_logEb1) {
        factor = m_a0*logE + m_b0;
    }
    else if (logE < m_logEb2) {
        factor = m_a1*logE + m_b1;
    }
    else {
        factor = m_a2*logE + m_b2;
    }

    // Return factor
    return factor;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear efficiency factor functor
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATEfficiency::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone efficiency factor functor
 *
 * @return Pointer to deep copy of efficiency factor functor.
 ***************************************************************************/
GLATEfficiency* GLATEfficiency::clone(void) const
{
    return new GLATEfficiency(*this);
}


/***********************************************************************//**
 * @brief Print efficiency factors
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing efficiency factors information.
 ***************************************************************************/
std::string GLATEfficiency::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATEfficiency ===");

        // Append information
        result.append("\n"+gammalib::parformat("Scale 1 (a0)")+gammalib::str(m_a0));
        result.append("\n"+gammalib::parformat("Scale 2 (a1)")+gammalib::str(m_a1));
        result.append("\n"+gammalib::parformat("Scale 3 (a2)")+gammalib::str(m_a2));
        result.append("\n"+gammalib::parformat("Offset 1 (b0)")+gammalib::str(m_b0));
        result.append("\n"+gammalib::parformat("Offset 2 (b1)")+gammalib::str(m_b1));
        result.append("\n"+gammalib::parformat("Offset 3 (b2)")+gammalib::str(m_b2));
        result.append("\n"+gammalib::parformat("Energy domains 1/2 limit (logEb1)")+gammalib::str(m_logEb1));
        result.append("\n"+gammalib::parformat("Energy domains 2/3 limit (logEb2)")+gammalib::str(m_logEb2));

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
void GLATEfficiency::init_members(void)
{
    // Initialise members
    m_a0     = 0.0;
    m_a1     = 0.0;
    m_a2     = 0.0;
    m_b0     = 0.0;
    m_b1     = 0.0;
    m_b2     = 0.0;
    m_logEb1 = 0.0;
    m_logEb2 = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] eff Efficiency factor functor.
 ***************************************************************************/
void GLATEfficiency::copy_members(const GLATEfficiency& eff)
{
    // Copy members
    m_a0     = eff.m_a0;
    m_a1     = eff.m_a1;
    m_a2     = eff.m_a2;
    m_b0     = eff.m_b0;
    m_b1     = eff.m_b1;
    m_b2     = eff.m_b2;
    m_logEb1 = eff.m_logEb1;
    m_logEb2 = eff.m_logEb2;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEfficiency::free_members(void)
{
    // Return
    return;
}
