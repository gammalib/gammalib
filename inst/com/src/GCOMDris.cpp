/***************************************************************************
 *             GCOMDris.cpp - COMPTEL Data Space container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2023 by Juergen Knodlseder                               *
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
 * @file GCOMDris.hpp
 * @brief COMPTEL Data Space container class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GMath.hpp"
#include "GCOMTools.hpp"
#include "GCOMSupport.hpp"
#include "GCOMDris.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMDris::at(int&)"
#define G_INSERT                           "GCOMDris::insert(int&, GCOMDri&)"
#define G_REMOVE                                     "GCOMDris::remove(int&)"

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
 *
 * Constructs an empty Data Space container.
 ***************************************************************************/
GCOMDris::GCOMDris(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dris Data Space container.
 ***************************************************************************/
GCOMDris::GCOMDris(const GCOMDris& dris)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(dris);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMDris::~GCOMDris(void)
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
 * @param[in] dris Data Space container.
 * @return Data Space container.
 ***************************************************************************/
GCOMDris& GCOMDris::operator=(const GCOMDris& dris)
{
    // Execute only if object is not identical
    if (this != &dris) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(dris);

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
 * @brief Clear Data Space container
 ***************************************************************************/
void GCOMDris::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone Data Space container
 *
 * @return Pointer to deep copy of Data Space container.
 ***************************************************************************/
GCOMDris* GCOMDris::clone(void) const
{
    return new GCOMDris(*this);
}


/***********************************************************************//**
 * @brief Return reference to Data Space
 *
 * @param[in] index Data Space index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Returns a reference to the Data Space with the specified @p index.
 ***************************************************************************/
GCOMDri& GCOMDris::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Data Space index",
                                       index, size());
    }

    // Return reference
    return m_dris[index];
}


/***********************************************************************//**
 * @brief Return reference to Data Space (const version)
 *
 * @param[in] index Data Space index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Returns a reference to the Data Space with the specified @p index.
 ***************************************************************************/
const GCOMDri& GCOMDris::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Data Space index",
                                       index, size());
    }

    // Return reference
    return m_dris[index];
}


/***********************************************************************//**
 * @brief Append Data Space to container
 *
 * @param[in] dri Data Space.
 * @return Reference to appended Data Space.
 *
 * Appends Data Space to the container by making a deep copy of the Data
 * Space.
 ***************************************************************************/
GCOMDri& GCOMDris::append(const GCOMDri& dri)
{
    // Append dri to list
    m_dris.push_back(dri);

    // Return reference
    return m_dris[size()-1];
}


/***********************************************************************//**
 * @brief Insert Data Space into container
 *
 * @param[in] index Data Space index (0,...,size()-1).
 * @param[in] dri Data Space.
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Inserts a Data Space into the container before the Data Space with the
 * specified @p index.
 ***************************************************************************/
GCOMDri& GCOMDris::insert(const int& index, const GCOMDri& dri)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Data Space index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Data Space index",
                                           index, size());
        }
    }
    #endif

    // Inserts Data Space
    m_dris.insert(m_dris.begin()+index, dri);

    // Return reference
    return m_dris[index];
}


/***********************************************************************//**
 * @brief Remove Data Space from container
 *
 * @param[in] index Data Space index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Data Space index is out of range.
 *
 * Remove Data Space of specified @p index from container.
 ***************************************************************************/
void GCOMDris::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Data Space index",
                                       index, size());
    }
    #endif

    // Erase Data Space from container
    m_dris.erase(m_dris.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Data Space container
 *
 * @param[in] oads Data Space container.
 *
 * Append Data Space container to the container.
 ***************************************************************************/
void GCOMDris::extend(const GCOMDris& dris)
{
    // Do nothing if Data Space container is empty
    if (!dris.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = dris.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_dris.push_back(dris[i]);
        }

    } // endif: Data Space container was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print Data Space container
 *
 * @param[in] chatter Chattiness.
 * @return String containing Data Space container information.
 ***************************************************************************/
std::string GCOMDris::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMDris ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Number of DRIs"));
        result.append(gammalib::str(size()));

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
void GCOMDris::init_members(void)
{
    // Initialise members
    m_dris.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dris Data Space container.
 ***************************************************************************/
void GCOMDris::copy_members(const GCOMDris& dris)
{
    // Copy members
    m_dris = dris.m_dris;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMDris::free_members(void)
{
    // Return
    return;
}
