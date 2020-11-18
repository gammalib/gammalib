/***************************************************************************
 *               GSkyDirs.cpp - Sky directions container class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSkyDirs.cpp
 * @brief Sky directions container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GSkyDirs.hpp"
/*
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableDoubleCol.hpp"
*/

/* __ Method name definitions ____________________________________________ */
#define G_AT                                            "GSkyDirs::at(int&)"
#define G_INSERT                          "GSkyDirs::insert(int&, GSkyDir&)"
#define G_REMOVE                                    "GSkyDirs::remove(int&)"

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
 * Constructs empty sky directions container.
 ***************************************************************************/
GSkyDirs::GSkyDirs(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Single sky direction constructor
 *
 * @param[in] dir Sky directions.
 *
 * Constructs sky directions container from a single sky direction.
 ***************************************************************************/
GSkyDirs::GSkyDirs(const GSkyDir& dir)
{
    // Initialise members
    init_members();

    // Append sky direction
    append(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] dir Sky directions.
 *
 * Construct sky directions container by copying from another sky directions
 * container.
 ***************************************************************************/
GSkyDirs::GSkyDirs(const GSkyDirs& dirs)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(dirs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GSkyDirs::~GSkyDirs(void)
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
 * @param[in] dirs Sky direction container.
 * @return Sky direction container.
 ***************************************************************************/
GSkyDirs& GSkyDirs::operator=(const GSkyDirs& dirs)
{
    // Execute only if object is not identical
    if (this != &dirs) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(dirs);

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
 * @brief Clear sky directions container
 *
 * Removes all sky directions from the container.
 ***************************************************************************/
void GSkyDirs::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone sky directions container
 *
 * @return Pointer to deep copy of sky directions container
 *
 * Makes a deep copy of the sky directions container instance.
 ***************************************************************************/
GSkyDirs* GSkyDirs::clone(void) const
{
    return new GSkyDirs(*this);
}


/***********************************************************************//**
 * @brief Return reference to sky direction (const version)
 *
 * @param[in] index Sky direction index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Sky direction index is out of range.
 *
 * Returns a reference to the sky direction with the specified @p index.
 ***************************************************************************/
const GSkyDir& GSkyDirs::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Sky direction", index, size());
    }

    // Return reference
    return m_dirs[index];
}


/***********************************************************************//**
 * @brief Append sky direction to container
 *
 * @param[in] dir Sky direction.
 * @return Reference to appended sky direction.
 *
 * Appends sky direction to the container by making a deep copy of the sky
 * direction.
 ***************************************************************************/
GSkyDir& GSkyDirs::append(const GSkyDir& dir)
{
    // Append sky direction to list
    m_dirs.push_back(dir);

    // Return reference
    return m_dirs[size()-1];
}


/***********************************************************************//**
 * @brief Insert sky direction into container
 *
 * @param[in] index Sky direction index (0,...,size()-1).
 * @param[in] dir Sky direction.
 *
 * @exception GException::out_of_range
 *            Sky direction index is out of range.
 *
 * Inserts a sky direction @p dir into the container before the sky
 * direction with the specified @p index.
 ***************************************************************************/
GSkyDir& GSkyDirs::insert(const int& index, const GSkyDir& dir)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "Sky direction", index,
                                           size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "Sky direction", index,
                                           size());
        }
    }
    #endif

    // Inserts sky direction
    m_dirs.insert(m_dirs.begin()+index, dir);

    // Return reference
    return m_dirs[index];
}


/***********************************************************************//**
 * @brief Remove sky direction from container
 *
 * @param[in] index Sky direction index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            Sky direction index is out of range.
 *
 * Remove sky direction of specified @p index from container.
 ***************************************************************************/
void GSkyDirs::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "Sky direction", index, size());
    }
    #endif

    // Erase sky direction from container
    m_dirs.erase(m_dirs.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append sky directions container
 *
 * @param[in] dirs Sky directions container.
 *
 * Append sky directions container to the container.
 ***************************************************************************/
void GSkyDirs::extend(const GSkyDirs& dirs)
{
    // Do nothing if sky directions container is empty
    if (!dirs.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = dirs.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_dirs.push_back(dirs[i]);
        }

    } // endif: sky directions container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print sky directions container information
 *
 * @param[in] chatter Chattiness.
 * @return String containing sky directions container information.
 ***************************************************************************/
std::string GSkyDirs::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GSkyDirs ===");

        // Append sky direction container information
        result.append("\n"+gammalib::parformat("Number of directions"));
        result.append(gammalib::str(size()));

        // EXPLICIT: Append sky directions
        if (chatter >= EXPLICIT) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n");
                result.append(gammalib::parformat("Direction "+gammalib::str(i)));
                result.append(m_dirs[i].print(chatter));
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
void GSkyDirs::init_members(void)
{
    // Initialise members
    m_dirs.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] dirs Sky directions container.
 ***************************************************************************/
void GSkyDirs::copy_members(const GSkyDirs& dirs)
{
    // Copy attributes
    m_dirs = dirs.m_dirs;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSkyDirs::free_members(void)
{
    // Return
    return;
}
