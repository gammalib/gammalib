/***************************************************************************
 *                  GEnergies.cpp - Energy container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GEnergies.cpp
 * @brief Energy container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GEnergies.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                            "GEnergies::at(int&)"


#define G_OP_ACCESS                                "GEnergies::operator[](int&)"
#define G_INSERT                               "GEnergies::insert(int&, GEnergy&)"
#define G_REMOVE                                       "GEnergies::remove(int&)"

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
GEnergies::GEnergies(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param energies Energy container.
 ***************************************************************************/
GEnergies::GEnergies(const GEnergies& energies)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(energies);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEnergies::~GEnergies(void)
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
 * @param[in] energies Energy container.
 * @return Energy container.
 ***************************************************************************/
GEnergies& GEnergies::operator=(const GEnergies& energies)
{
    // Execute only if object is not identical
    if (this != &energies) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(energies);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear container
 *
 * Removes all energies from the container.
 ***************************************************************************/
void GEnergies::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of energy container
 *
 * Makes a deep copy of the energy container instance.
 ***************************************************************************/
GEnergies* GEnergies::clone(void) const
{
    return new GEnergies(*this);
}


/***********************************************************************//**
 * @brief Return reference to energy
 *
 * @param[in] index Energy index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Returns a reference to the energy with the specified @p index.
 ***************************************************************************/
GEnergy& GEnergies::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return m_energies[index];
}


/***********************************************************************//**
 * @brief Return reference to energy (const version)
 *
 * @param[in] index Energy index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Returns a reference to the energy with the specified @p index.
 ***************************************************************************/
const GEnergy& GEnergies::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, index, 0, size()-1);
    }

    // Return reference
    return m_energies[index];
}


/***********************************************************************//**
 * @brief Append energy to container
 *
 * @param[in] energy Energy.
 * @return Reference to appended energy.
 *
 * Appends energy to the container by making a deep copy of the energy.
 ***************************************************************************/
GEnergy& GEnergies::append(const GEnergy& energy)
{
    // Append energy to list
    m_energies.push_back(energy);

    // Return reference
    return m_energies[size()-1];
}


/***********************************************************************//**
 * @brief Insert energy into container
 *
 * @param[in] index Energy index (0,...,size()-1).
 * @param[in] energy Energy.
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Inserts an @p energy into the container before the energy with the
 * specified @p index.
 ***************************************************************************/
GEnergy& GEnergies::insert(const int& index, const GEnergy& energy)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (isempty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, index, 0, size()-1);
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, index, 0, size()-1);
        }
    }
    #endif

    // Inserts energy
    m_energies.insert(m_energies.begin()+index, energy);

    // Return reference
    return m_energies[index];
}


/***********************************************************************//**
 * @brief Remove energy from container
 *
 * @param[in] index Energy index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Energy index is out of range.
 *
 * Remove energy of specified @p index from container.
 ***************************************************************************/
void GEnergies::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, index, 0, size()-1);
    }
    #endif

    // Erase energy from container
    m_energies.erase(m_energies.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append energy container
 *
 * @param[in] energies Energy container.
 *
 * Append energy container to the container.
 ***************************************************************************/
void GEnergies::extend(const GEnergies& energies)
{
    // Do nothing if energy container is empty
    if (!energies.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = energies.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_energies.push_back(energies[i]);
        }

    } // endif: energy container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print energy container information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy container information.
 ***************************************************************************/
std::string GEnergies::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GEnergies ===");

        // Append energy container information
        result.append("\n"+gammalib::parformat("Number of energies"));
        result.append(gammalib::str(size()));

        // EXPLICIT: Append energies
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+m_energies[i].print(chatter));
            }
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEnergies::init_members(void)
{
    // Initialise members
    m_energies.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] energies Energy container.
 ***************************************************************************/
void GEnergies::copy_members(const GEnergies& energies)
{
    // Copy attributes
    m_energies = energies.m_energies;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEnergies::free_members(void)
{
    // Return
    return;
}
