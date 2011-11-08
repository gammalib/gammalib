/***************************************************************************
 *         GObservationRegistry.cpp  -  Observation registry class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Juergen Knoedlseder                              *
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
 * @brief GObservationRegistry class interface implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservationRegistry.hpp"
#include "GTools.hpp"

/* __ Static members _____________________________________________________ */
int                  GObservationRegistry::m_number(0);
std::string*         GObservationRegistry::m_names(0);
const GObservation** GObservationRegistry::m_obs(0);

/* __ Method name definitions ____________________________________________ */
#define G_INSTRUMENT                 "GObservationRegistry::instrument(int&)"

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
    for (int i = 0; i < m_number; ++i) {
        std::cout << "\"" << m_names[i] << "\" ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Observation constructor
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
    std::string*         new_names = new std::string[m_number+1];
    const GObservation** new_obs   = new const GObservation*[m_number+1];

    // Save old registry
    for (int i = 0; i < m_number; ++i) {
        new_names[i] = m_names[i];
        new_obs[i]   = m_obs[i];
    }

    // Add new observation to registry
    new_names[m_number] = obs->instrument();
    new_obs[m_number]   = obs;

    // Delete old registry
    if (m_names != NULL) delete [] m_names;
    if (m_obs   != NULL) delete [] m_obs;

    // Set pointers on new registry
    m_names = new_names;
    m_obs   = new_obs;

    // Increment number of observations in registry
    m_number++;

    // Debug option: Show actual registry
    #if G_DEBUG_REGISTRY
    std::cout << "GObservationRegistry(const GObservation*): ";
    for (int i = 0; i < m_number; ++i) {
        std::cout << "\"" << m_names[i] << "\" ";
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
 * @brief Allocate observation of given type
 *
 * @param[in] instrument Instrument name.
 *
 * Returns a pointer to a void observation instance for a specific
 * instrument. If the instrument has not been found in the registry, a NULL
 * pointer is returned.
 ***************************************************************************/
GObservation* GObservationRegistry::alloc(const std::string& instrument) const
{
    // Initialise observation
    GObservation* obs = NULL;

    // Search for observation in registry
    for (int i = 0; i < m_number; ++i) {
        if (m_names[i] == instrument) {
            obs = m_obs[i]->clone();
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
 *
 * @exception GException::out_of_range
 *            Observation index is out of range.
 ***************************************************************************/
std::string GObservationRegistry::instrument(const int& index) const
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_INSTRUMENT, index, 0, size()-1);
    }
    #endif

    // Return name
    return (m_names[index]);
}


/***********************************************************************//**
 * @brief Print registry information
 ***************************************************************************/
std::string GObservationRegistry::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GObservationRegistry ===");
    result.append("\n"+parformat("Number of observations")+str(m_number));

    // Append observations
    for (int i = 0; i < m_number; ++i) {
        result.append("\n"+parformat(m_names[i]));
        result.append(m_obs[i]->instrument());
    }

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


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] registry Observation registry.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GObservationRegistry& registry)
{
     // Write registry in output stream
    os << registry.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] registry Observation registry.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GObservationRegistry& registry)
{
    // Write registry into logger
    log << registry.print();

    // Return logger
    return log;
}
