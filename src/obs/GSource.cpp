/***************************************************************************
 *                        GSource.cpp - Source class                       *
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
 * @file GSource.cpp
 * @brief Source class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSource.hpp"
#include "GTools.hpp"

/* __ Constants __________________________________________________________ */

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
GSource::GSource(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Source constructor
 *
 * @param[in] name Source name.
 * @param[in] model Spatial model pointer.
 * @param[in] energy Energy.
 * @param[in] time Time.
 ***************************************************************************/
GSource::GSource(const std::string& name,
                 GModelSpatial*     model,
                 const GEnergy&     energy,
                 const GTime&       time)
{ 
    // Initialise private members
    init_members();

    // Set members
    m_name   = name;
    m_model  = model;
    m_energy = energy;
    m_time   = time;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] src Source.
 ***************************************************************************/
GSource::GSource(const GSource& src)
{ 
    // Initialise private members
    init_members();

    // Copy members
    copy_members(src);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSource::~GSource(void)
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
 * @param[in] src Source.
 * @return Source.
 ***************************************************************************/
GSource& GSource::operator= (const GSource& src)
{ 
    // Execute only if object is not identical
    if (this != &src) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(src);

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
 * @brief Clear instance
 ***************************************************************************/
void GSource::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pointer to deep copy of source.
 ***************************************************************************/
GSource* GSource::clone(void) const
{
    // Clone this image
    return new GSource(*this);
}


/***********************************************************************//**
 * @brief Print source
 *
 * @return String containing source information
 ***************************************************************************/
std::string GSource::print(void) const
{
    // Initialise result string
    std::string result;

    // Build photon string
    result.append("GSource(");
    result.append(m_name);
    result.append(", "+m_model->type());
    result.append(", E="+m_energy.print());
    result.append(", MET="+m_time.print());
    result.append(")");

    // Return
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
void GSource::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_model = NULL;
    m_energy.clear();
    m_time.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] src Source.
 ***************************************************************************/
void GSource::copy_members(const GSource& src)
{
    // Copy members
    m_model  = src.m_model;
    m_name   = src.m_name;
    m_energy = src.m_energy;
    m_time   = src.m_time;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSource::free_members(void)
{
    // Return
    return;
}
