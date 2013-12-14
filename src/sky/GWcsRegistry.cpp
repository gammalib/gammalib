/***************************************************************************
 *      GWcsRegistry.cpp - World Coordinate Projection registry class      *
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
 * @file GWcsRegistry.cpp
 * @brief World Coordinate Projection registry class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GWcsRegistry.hpp"
#include "GException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CODE                                     "GWcsRegistry::code(int&)"
#define G_NAME                                     "GWcsRegistry::name(int&)"

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
GWcsRegistry::GWcsRegistry(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GWcsRegistry(void): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << codes()[i] << "\" ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Projection constructor
 *
 * @param[in] wcs World Coordinate System.
 ***************************************************************************/
GWcsRegistry::GWcsRegistry(const GWcs* wcs)
{
    // Initialise private members for clean destruction
    init_members();

    // Debug option: Notify new registry
    #if G_DEBUG_REGISTRY
    std::cout << "GWcsRegistry(const GWcs*): ";
    std::cout << "add \"" << wcs->code() << "\" to registry." << std::endl;
    #endif

    // Allocate new registry
    std::string* new_codes        = new std::string[size()+1];
    std::string* new_names        = new std::string[size()+1];
    const GWcs** new_projections  = new const GWcs*[size()+1];

    // Save old registry
    for (int i = 0; i < size(); ++i) {
        new_codes[i]       = codes()[i];
        new_names[i]       = names()[i];
        new_projections[i] = projections()[i];
    }

    // Add new projection to registry
    new_codes[size()]       = wcs->code();
    new_names[size()]       = wcs->name();
    new_projections[size()] = wcs;

    // Set pointers on new registry
    codes().assign(new_codes);
    names().assign(new_names);
    projections().assign(new_projections);

    // Increment number of projections in registry
    number()++;

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GWcsRegistry(const GWcs*): ";
    for (int i = 0; i < size(); ++i) {
        std::cout << "\"" << codes()[i] << "\" ";
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
GWcsRegistry::GWcsRegistry(const GWcsRegistry& registry)
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
GWcsRegistry::~GWcsRegistry(void)
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
GWcsRegistry& GWcsRegistry::operator=(const GWcsRegistry& registry)
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
 * @brief Allocate World Coordinate System of given code
 *
 * @param[in] code Sky projection code.
 * @return Pointer to World Coordinate System (NULL if code is not registered).
 *
 * Returns a pointer to a sky projection instance of the specified code. If
 * the code has not been found in the registry, a NULL pointer is returned.
 ***************************************************************************/
GWcs* GWcsRegistry::alloc(const std::string& code) const
{
    // Initialise projection
    GWcs* projection = NULL;

    // Search for projection in registry
    for (int i = 0; i < size(); ++i) {
        if (codes()[i] == code) {
            projection = projections()[i]->clone();
            break;
        }
    }

    // Return projection
    return projection;
}


/***********************************************************************//**
 * @brief Returns projection code
 *
 * @param[in] index Projection index [0,...,size()-1].
 * @return Projection code.
 *
 * @exception GException::out_of_range
 *            Projection index is out of range.
 ***************************************************************************/
std::string GWcsRegistry::code(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_CODE, index, 0, size()-1);
    }
    #endif

    // Return code
    return (codes()[index]);
}


/***********************************************************************//**
 * @brief Returns projection name
 *
 * @param[in] index Projection index [0,...,size()-1].
 * @return Projection name.
 *
 * @exception GException::out_of_range
 *            Projection index is out of range.
 ***************************************************************************/
std::string GWcsRegistry::name(const int& index) const
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
 * @brief Return list string of projection codes
 *
 * @return Projection code list.
 *
 * Returns a list of the projection codes in the format 'xxx/yyy/zzz'.
 ***************************************************************************/
std::string GWcsRegistry::list(void) const
{
    // Initialise result string
    std::string result;

    // Append projections
    for (int i = 0; i < size(); ++i) {
        if (i > 0) {
            result.append("/");
        }
        result.append(codes()[i]);
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print registry information
 *
 * @return Registry content.
 ***************************************************************************/
std::string GWcsRegistry::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GWcsRegistry ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of projections"));
        result.append(gammalib::str(size()));

        // NORMAL: Append projections
        if (chatter >= NORMAL) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+gammalib::parformat(codes()[i]));
                result.append(names()[i]);
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
void GWcsRegistry::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] registry Registry.
 ***************************************************************************/
void GWcsRegistry::copy_members(const GWcsRegistry& registry)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GWcsRegistry::free_members(void)
{
    // Return
    return;
}
