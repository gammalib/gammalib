/***************************************************************************
 *                GTPLContainer.cpp - [WHAT] container class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GTPLContainer.hpp
 * @brief [WHAT] container class implementation
 * @author [AUTHOR]
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTPLContainer.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_AT                                        "GTPLContainer::at(int&)"
#define G_INSERT                     "GTPLContainer::insert(int&, GTPLBase&)"
#define G_REMOVE                                "GTPLContainer::remove(int&)"

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
GTPLContainer::GTPLContainer(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] TPL_CONTAINER [WHAT] container.
 ***************************************************************************/
GTPLContainer::GTPLContainer(const GTPLContainer& TPL_CONTAINER)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(TPL_CONTAINER);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTPLContainer::~GTPLContainer(void)
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
 * @param[in] TPL_CONTAINER [WHAT] container.
 * @return [WHAT] container.
 ***************************************************************************/
GTPLContainer& GTPLContainer::operator=(const GTPLContainer& TPL_CONTAINER)
{
    // Execute only if object is not identical
    if (this != &TPL_CONTAINER) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(TPL_CONTAINER);

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
 * @brief Clear [WHAT] container
 ***************************************************************************/
void GTPLContainer::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone [WHAT] container
 *
 * @return Pointer to deep copy of [WHAT] container.
 ***************************************************************************/
GTPLContainer* GTPLContainer::clone(void) const
{
    return new GTPLContainer(*this);
}


/***********************************************************************//**
 * @brief Return reference to [WHAT]
 *
 * @param[in] index [WHAT] index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            [WHAT] index is out of range.
 *
 * Returns a reference to the [WHAT] with the specified @p index.
 ***************************************************************************/
GTPLBase& GTPLContainer::at(const int& index)
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "[WHAT] index",
                                       index, size());
    }

    // Return reference
    return m_TPL_CONTAINER[index];
}


/***********************************************************************//**
 * @brief Return reference to [WHAT] (const version)
 *
 * @param[in] index [WHAT] index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            [WHAT] index is out of range.
 *
 * Returns a reference to the [WHAT] with the specified @p index.
 ***************************************************************************/
const GTPLBase& GTPLContainer::at(const int& index) const
{
    // Raise exception if index is out of range
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "[WHAT] index",
                                       index, size());
    }

    // Return reference
    return m_TPL_CONTAINER[index];
}


/***********************************************************************//**
 * @brief Append [WHAT] to container
 *
 * @param[in] TPL_OBJECT [WHAT].
 * @return Reference to appended [WHAT].
 *
 * Appends [WHAT] to the container by making a deep copy of the
 * [WHAT].
 ***************************************************************************/
GTPLBase& GTPLContainer::append(const GTPLBase& TPL_OBJECT)
{
    // Append TPL_OBJECT to list
    m_TPL_CONTAINER.push_back(TPL_OBJECT);

    // Return reference
    return m_TPL_CONTAINER[size()-1];
}


/***********************************************************************//**
 * @brief Insert [WHAT] into container
 *
 * @param[in] index [WHAT] index (0,...,size()-1).
 * @param[in] TPL_OBJECT [WHAT].
 *
 * @exception GException::out_of_range
 *            [WHAT] index is out of range.
 *
 * Inserts an @p TPL_OBJECT into the container before the [WHAT] with the
 * specified @p index.
 ***************************************************************************/
GTPLBase& GTPLContainer::insert(const int& index, const GTPLBase& TPL_OBJECT)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT, "[WHAT] index",
                                           index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT, "[WHAT] index",
                                           index, size());
        }
    }
    #endif

    // Inserts [WHAT]
    m_TPL_CONTAINER.insert(m_TPL_CONTAINER.begin()+index, TPL_OBJECT);

    // Return reference
    return m_TPL_CONTAINER[index];
}


/***********************************************************************//**
 * @brief Remove [WHAT] from container
 *
 * @param[in] index [WHAT] index (0,...,size()-1).
 *
 * @exception GException::out_of_range
 *            [WHAT] index is out of range.
 *
 * Remove [WHAT] of specified @p index from container.
 ***************************************************************************/
void GTPLContainer::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE, "[WHAT] index",
                                       index, size());
    }
    #endif

    // Erase [WHAT] from container
    m_TPL_CONTAINER.erase(m_TPL_CONTAINER.begin() + index);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Append [WHAT] container
 *
 * @param[in] TPL_CONTAINER [WHAT] container.
 *
 * Append [WHAT] container to the container.
 ***************************************************************************/
void GTPLContainer::extend(const GTPLContainer& TPL_CONTAINER)
{
    // Do nothing if [WHAT] container is empty
    if (!TPL_CONTAINER.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = TPL_CONTAINER.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all elements and append them to container
        for (int i = 0; i < num; ++i) {
            m_TPL_CONTAINER.push_back(TPL_CONTAINER[i]);
        }

    } // endif: [WHAT] container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print [WHAT] container
 *
 * @param[in] chatter Chattiness.
 * @return String containing [WHAT] container information.
 *
 * @todo Implement method.
 ***************************************************************************/
std::string GTPLContainer::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GTPLContainer ===");

        // Append container information
        result.append("\n"+gammalib::parformat("Number of TPL_OBJECT"));
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
void GTPLContainer::init_members(void)
{
    // Initialise members
    m_TPL_CONTAINER.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] TPL_OBJECT [WHAT] container.
 ***************************************************************************/
void GTPLContainer::copy_members(const GTPLContainer& TPL_CONTAINER)
{
    // Copy members
    m_TPL_CONTAINER = TPL_CONTAINER.m_TPL_CONTAINER;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTPLContainer::free_members(void)
{
    // Return
    return;
}
