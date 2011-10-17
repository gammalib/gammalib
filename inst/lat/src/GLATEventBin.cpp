/***************************************************************************
 *                GLATEventBin.cpp  -  LAT event bin class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
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
 * @file GLATEventBin.cpp
 * @brief GLATEventBin class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <cmath>
#include "GLATEventBin.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATEventBin::GLATEventBin(void) : GEventBin()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
GLATEventBin::GLATEventBin(const GLATEventBin& bin) : GEventBin(bin)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(bin);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATEventBin::~GLATEventBin(void)
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
 * @param[in] bin Event bin which should be assigned.
 ***************************************************************************/
GLATEventBin& GLATEventBin::operator= (const GLATEventBin& bin)
{
    // Execute only if object is not identical
    if (this != &bin) {

        // Copy base class members
        this->GEventBin::operator=(bin);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(bin);

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
 * @brief Clear instance
 *
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GLATEventBin::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GEventBin::free_members();
    this->GEvent::free_members();

    // Initialise members
    this->GEvent::init_members();
    this->GEventBin::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATEventBin* GLATEventBin::clone(void) const
{
    return new GLATEventBin(*this);
}


/***********************************************************************//**
 * @brief Return size of event bin
 *
 * The size of the event bin (units sr MeV s) is given by
 * \f[size = \Omega \times \Delta E \times \Delta T\f]
 * where
 * \f$\Omega\f$ is the size of the spatial bin in sr,
 * \f$\Delta E\f$ is the size of the energy bin in MeV, and
 * \f$\Delta T\f$ is the ontime of the observation in seconds. 
 ***************************************************************************/
double GLATEventBin::size(void) const
{
    // Compute bin size
    double size = omega() * ewidth().MeV() * ontime();

    // Return bin size
    return size;
}


/***********************************************************************//**
 * @brief Return error in number of counts
 *
 * Returns \f$\sqrt(counts+delta)\f$ as the uncertainty in the number of
 * counts in the bin. Adding delta avoids uncertainties of 0 which will
 * lead in the optimisation step to the exlusion of the corresponding bin.
 * In the actual implementation delta=1e-50.
 *
 * @todo The choice of delta has been made somewhat arbitrary, mainly
 * because the optimizer routines filter error^2 below 1e-100.
 ***************************************************************************/
double GLATEventBin::error(void) const
{
    // Compute uncertainty
    double error = sqrt(counts()+1.0e-50);

    // Return error
    return error;
}


/***********************************************************************//**
 * @brief Print event information
 ***************************************************************************/
std::string GLATEventBin::print(void) const
{
    // Initialise result string
    std::string result;

    // Append number of counts
    result.append(str(counts()));

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATEventBin::init_members(void)
{
    // Initialise members
    m_cube   = NULL;
    m_index  = 0;
    m_ipix   = 0;
    m_ieng   = 0;
    m_energy = NULL;
    m_dir    = NULL;
    m_time   = NULL;
    m_counts = NULL;
    m_omega  = NULL;
    m_ewidth = NULL;
    m_ontime = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] bin Event bin.
 ***************************************************************************/
void GLATEventBin::copy_members(const GLATEventBin& bin)
{
    // Copy members
    m_cube   = bin.m_cube;
    m_index  = bin.m_index;
    m_ipix   = bin.m_ipix;
    m_ieng   = bin.m_ieng;
    m_energy = bin.m_energy;
    m_dir    = bin.m_dir;
    m_time   = bin.m_time;
    m_counts = bin.m_counts;
    m_omega  = bin.m_omega;
    m_ewidth = bin.m_ewidth;
    m_ontime = bin.m_ontime;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATEventBin::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
