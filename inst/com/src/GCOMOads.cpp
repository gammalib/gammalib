/***************************************************************************
 *         GCOMOads.cpp - COMPTEL Orbit Aspect Data container class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMOads.hpp
 * @brief COMPTEL Orbit Aspect Data container class implementation
 * @author Juergen Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GCOMOads.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                             "GCOMOads::at(int&)"
#define G_INSERT                           "GCOMOads::insert(int&, GCOMOad&)"
#define G_REMOVE                                     "GCOMOads::remove(int&)"

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
GCOMOads::GCOMOads(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
GCOMOads::GCOMOads(const GCOMOads& oads)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(oads);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCOMOads::~GCOMOads(void)
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
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 * @return COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
GCOMOads& GCOMOads::operator=(const GCOMOads& oads)
{
    // Execute only if object is not identical
    if (this != &oads) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(oads);

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
 * @brief Clear COMPTEL Orbit Aspect Data container
 ***************************************************************************/
void GCOMOads::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone COMPTEL Orbit Aspect Data container
 *
 * @return Pointer to deep copy of COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
GCOMOads* GCOMOads::clone(void) const
{
    return new GCOMOads(*this);
}


/***********************************************************************//**
 * @brief Return reference to Orbit Aspect Data
 *
 * @param[in] index Orbit Aspect Data index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Returns a reference to the Orbit Aspect Data with the specified @p index.
 ***************************************************************************/
GCOMOad& GCOMOads::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Orbit Aspect Data index",
                                       index, size());
    }

    // Return reference
    return m_oads[index];
}


/***********************************************************************//**
 * @brief Return reference to Orbit Aspect Data (const version)
 *
 * @param[in] index Orbit Aspect Data index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Returns a reference to the Orbit Aspect Data with the specified @p index.
 ***************************************************************************/
const GCOMOad& GCOMOads::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Orbit Aspect Data index",
                                       index, size());
    }

    // Return reference
    return m_oads[index];
}


/***********************************************************************//**
 * @brief Append Orbit Aspect Data to container
 *
 * @param[in] oad Orbit Aspect Data.
 * @return Reference to appended Orbit Aspect Data.
 *
 * Appends Orbit Aspect Data to the container by making a deep copy of the
 * Orbit Aspect Data.
 ***************************************************************************/
GCOMOad& GCOMOads::append(const GCOMOad& oad)
{
    // Append oad to list
    m_oads.push_back(oad);

    // Return reference
    return m_oads[size()-1];
}


/***********************************************************************//**
 * @brief Insert Orbit Aspect Data into container
 *
 * @param[in] index Orbit Aspect Data index (0,...,size()-1).
 * @param[in] oad Orbit Aspect Data.
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Inserts an @p Orbit Aspect Data into the container before the Orbit
 * Aspect Data with the specified @p index.
 ***************************************************************************/
GCOMOad& GCOMOads::insert(const int& index, const GCOMOad& oad)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Orbit Aspect Data index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Orbit Aspect Data index",
                                           index, size());
        }
    }
    #endif

    // Inserts Orbit Aspect Data
    m_oads.insert(m_oads.begin()+index, oad);

    // Return reference
    return m_oads[index];
}


/***********************************************************************//**
 * @brief Remove Orbit Aspect Data from container
 *
 * @param[in] index Orbit Aspect Data index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Orbit Aspect Data index is out of range.
 *
 * Remove Orbit Aspect Data of specified @p index from container.
 ***************************************************************************/
void GCOMOads::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Orbit Aspect Data index",
                                       index, size());
    }
    #endif

    // Erase Orbit Aspect Data from container
    m_oads.erase(m_oads.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append Orbit Aspect Data container
 *
 * @param[in] oads COMPTEL Orbit Aspect Data container.
 *
 * Append COMPTEL Orbit Aspect Data container to the container.
 ***************************************************************************/
void GCOMOads::extend(const GCOMOads& oads)
{
    // Do nothing if COMPTEL Orbit Aspect Data container is empty
    if (!oads.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = oads.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_oads.push_back(oads[i]);
        }

    } // endif: COMPTEL Orbit Aspect Data container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL Orbit Aspect Data container
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL Orbit Aspect Data container information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GCOMOads::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMOads ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Number of oad"));
        result.append(gammalib::str(size()));

        // Append other information
        // TODO: Add any relevant information

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
void GCOMOads::init_members(void)
{
    // Initialise members
    m_oads.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] oad COMPTEL Orbit Aspect Data container.
 ***************************************************************************/
void GCOMOads::copy_members(const GCOMOads& oads)
{
    // Copy members
    m_oads = oads.m_oads;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMOads::free_members(void)
{
    // Return
    return;
}
