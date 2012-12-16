/***************************************************************************
 *                     GTimes.cpp - Time container class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GTimes.cpp
 * @brief Time container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GTimes.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS                                "GTimes::operator[](int&)"

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
GTimes::GTimes(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param times Photon container.
 ***************************************************************************/
GTimes::GTimes(const GTimes& times)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(times);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTimes::~GTimes(void)
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
 * @param[in] times Time container.
 * @return Time container.
 ***************************************************************************/
GTimes& GTimes::operator=(const GTimes& times)
{
    // Execute only if object is not identical
    if (this != &times) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(times);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to time
 *
 * @param[in] index Index of time [0,...,size()-1]
 *
 * @exception GException::out_of_range
 *            Time index is out of range.
 ***************************************************************************/
GTime& GTimes::operator[](const int& index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OP_ACCESS, index, 0, size()-1);
    }

    // Return reference
    return m_times[index];
}


/***********************************************************************//**
 * @brief Return reference to time (const version)
 *
 * @param[in] index Index of time [0,...,size()-1]
 *
 * @exception GException::out_of_range
 *            Time index is out of range.
 ***************************************************************************/
const GTime& GTimes::operator[](const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_OP_ACCESS, index, 0, size()-1);
    }

    // Return reference
    return m_times[index];
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear container
 ***************************************************************************/
void GTimes::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 *
 * @return Pointer to deep copy of time container.
 ***************************************************************************/
GTimes* GTimes::clone(void) const
{
    // Clone this image
    return new GTimes(*this);
}


/***********************************************************************//**
 * @brief Append time to container
 *
 * @param[in] time Time.
 *
 * This method appends a time to the container by copying it.
 ***************************************************************************/
void GTimes::append(const GTime& time)
{
    // Append time to list
    m_times.push_back(time);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve memory for times in container
 *
 * @param[in] number Number of times.
 *
 * This method reserves memory for a number of times in the container.
 ***************************************************************************/
void GTimes::reserve(const int& number)
{
    // Reserve memory
    m_times.reserve(number);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print time container information
 *
 * @return String containing time container information.
 ***************************************************************************/
std::string GTimes::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GTimes ===\n");

    // Append time container information
    result.append(parformat("Number of times")+str(size())+"\n");

    // Append times
    for (int i = 0; i < size(); ++i) {
        result.append("\n");
        result.append(m_times[i].print());
    }

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
void GTimes::init_members(void)
{
    // Initialise members
    m_times.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] times Time container.
 ***************************************************************************/
void GTimes::copy_members(const GTimes& times)
{
    // Copy attributes
    m_times = times.m_times;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTimes::free_members(void)
{
    // Return
    return;
}
