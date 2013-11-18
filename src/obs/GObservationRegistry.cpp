/***************************************************************************
 *          GObservationRegistry.cpp - Observation registry class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GObservationRegistry.cpp
 * @brief Observation registry class definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservationRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_NAME                             "GObservationRegistry::name(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_REGISTRY 0


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GObservationRegistry::GObservationRegistry(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GObservationRegistry(void): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << names()[i] << "\" ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observation constructor
 *
 * @param[in] obs Observation.
 ***************************************************************************/
GObservationRegistry::GObservationRegistry(const GObservation* obs)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GObservationRegistry(const GObservation*): ";
    std::cout << "add \"" << obs->instrument() << "\" to registry." << std::endl;
    #endif

    // Allocate new registry
    std::string*         new_names = new std::string[size()+1];
    const GObservation** new_obs   = new const GObservation*[size()+1];

    // Save old registry
    for (int i = 0; i < size(); ++i) {
        new_names[i] = names()[i];
        new_obs[i]   = this->obs()[i];
    }

    // Add new observation to registry
    new_names[size()] = obs->instrument();
    new_obs[size()]   = obs;

    // Set pointers on new registry
    names().assign(new_names);
    this->obs().assign(new_obs);

    // Increment number of observations in registry
    number()++;

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GObservationRegistry(const GObservation*): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << names()[i] << "\" ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] registry Registry.
 ***************************************************************************/
GObservationRegistry::GObservationRegistry(const GObservationRegistry& registry)
{
    // Initialise private members
    init_members();

    // Copy members
    copy_members(registry);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GObservationRegistry::~GObservationRegistry(void)
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
 * @param[in] registry Registry.
 * @return Reference to registry.
 ***************************************************************************/
GObservationRegistry& GObservationRegistry::operator= (const GObservationRegistry& registry)
{
    // Execute only if object is not identical
    if (this != &registry) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(registry);

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
 * @brief Allocate observation of given name
 *
 * @param[in] name Instrument name.
 * @return Pointer to observation (NULL if name is not registered).
 *
 * Returns a pointer to an observation instance for a specific name.
 * If the instrument name has not been found in the registry, a NULL
 * pointer is returned.
 ***************************************************************************/
GObservation* GObservationRegistry::alloc(const std::string& name) const
{
    // Initialise observation
    GObservation* obs = NULL;

    // Search for observation in registry
    for (int i = 0; i < size(); ++i) {
        if (names()[i] == name) {
            obs = this->obs()[i]->clone();
            break;
        }
    }

    // Return observation
    return obs;
}


/***********************************************************************//**
 * @brief Returns instrument name for a specific registered observation
 *
 * @param[in] index Observation index [0,...,size()-1].
 * @return Instrument name.
 *
 * @exception GException::out_of_range
 *            Observation index is out of range.
 ***************************************************************************/
std::string GObservationRegistry::name(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_NAME, index, 0, size()-1);
    }
    #endif

    // Return name
    return (names()[index]);
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return Registry content.
 ***************************************************************************/
std::string GObservationRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GObservationRegistry ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of observations"));
        result.append(gammalib::str(size()));

        // NORMAL: Append observations
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+gammalib::parformat(names()[i]));
                result.append(this->obs()[i]->instrument());
            }
        }

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
void GObservationRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GObservationRegistry::copy_members(const GObservationRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GObservationRegistry::free_members(void)
{
    // Return
    return;
}
