/***************************************************************************
 *                GTPLContainer.hpp - [WHAT] container class               *
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
 * @brief [WHAT] container class definition
 * @author [AUTHOR]
 */

#ifndef GTPLCONTAINER_HPP
#define GTPLCONTAINER_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GTPLBase.hpp"
#include "GContainer.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GTPLContainer
 *
 * @brief [WHAT] container class
 *
 * @todo Add class description.
 ***************************************************************************/
class GTPLContainer : public GContainer {

public:
    // Constructors and destructors
    GTPLContainer(void);
    GTPLContainer(const GTPLContainer& TPL_CONTAINER);
    virtual ~GTPLContainer(void);

    // Operators
    GTPLContainer&  operator=(const GTPLContainer& TPL_CONTAINER);
    GTPLBase&       operator[](const int& index);
    const GTPLBase& operator[](const int& index) const;

    // Methods
    void            clear(void);
    GTPLContainer*  clone(void) const;
    std::string     classname(void) const;
    GTPLBase&       at(const int& index);
    const GTPLBase& at(const int& index) const;
    int             size(void) const;
    bool            is_empty(void) const;
    GTPLBase&       append(const GTPLBase& TPL_OBJECT);
    GTPLBase&       insert(const int& index, const GTPLBase& TPL_OBJECT);
    void            remove(const int& index);
    void            reserve(const int& num);
    void            extend(const GTPLContainer& TPL_CONTAINER);
    std::string     print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTPLContainer& TPL_CONTAINER);
    void free_members(void);
    
    // Protected data members
    std::vector<GTPLBase> m_TPL_CONTAINER;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GTPLContainer").
 ***************************************************************************/
inline
std::string GTPLContainer::classname(void) const
{
    return ("GTPLContainer");
}


/***********************************************************************//**
 * @brief Return reference to TPL_OBJECT
 *
 * @param[in] index TPL_OBJECT index [0,...,size()-1].
 *
 * Returns a reference to the TPL_OBJECT with the specified @p index.
 ***************************************************************************/
inline
GTPLBase& GTPLContainer::operator[](const int& index)
{
    return (m_TPL_CONTAINER[index]);
}


/***********************************************************************//**
 * @brief Return reference to TPL_OBJECT (const version)
 *
 * @param[in] index TPL_OBJECT index [0,...,size()-1].
 *
 * Returns a reference to the TPL_OBJECT with the specified @p index.
 ***************************************************************************/
inline
const GTPLBase& GTPLContainer::operator[](const int& index) const
{
    return (m_TPL_CONTAINER[index]);
}


/***********************************************************************//**
 * @brief Return number of TPL_CONTAINER in container
 *
 * @return Number of TPL_CONTAINER in container.
 *
 * Returns the number of TPL_CONTAINER in the container.
 ***************************************************************************/
inline
int GTPLContainer::size(void) const
{
    return (int)m_TPL_CONTAINER.size();
}


/***********************************************************************//**
 * @brief Signals if there are no TPL_CONTAINER in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the TPL_OBJECT container does not contain any TPL_OBJECT.
 ***************************************************************************/
inline
bool GTPLContainer::is_empty(void) const
{
    return (m_TPL_CONTAINER.empty());
}


/***********************************************************************//**
 * @brief Reserves space for TPL_CONTAINER in container
 *
 * @param[in] num Number of TPL_CONTAINER.
 *
 * Reserves space for @p num TPL_CONTAINER in the container.
 ***************************************************************************/
inline
void GTPLContainer::reserve(const int& num)
{
    m_TPL_CONTAINER.reserve(num);
    return;
}

#endif /* GTPLCONTAINER_HPP */
