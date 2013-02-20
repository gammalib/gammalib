/***************************************************************************
 *                GContainer.hpp - Container interface class               *
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
 * @file GContainer.hpp
 * @brief Definition of interface for container classes
 * @author Juergen Knoedlseder
 */

#ifndef GCONTAINER_HPP
#define GCONTAINER_HPP

/* __ Includes ___________________________________________________________ */
#include <GBase.hpp>


/***********************************************************************//**
 * @class GContainer
 *
 * @brief Interface class for container classes
 *
 * This class defines the interface for container classes. The usage of the
 * interface class imposes on all container classes are coherent interface.
 * The following methods are mandatory:
 *
 *     clear      - Clear container (inherited from GBase)
 *     clone      - Clones container (inherited from GBase)
 *     print      - Print container content (inherited from GBase)
 *     size       - Returns number of objects is container
 *     isempty    - Checks if container is empty
 *     remove     - Removes an object from the container
 *     reserve    - Reserves space in the container
 ***************************************************************************/
class GContainer : public GBase {

public:
    /// @brief Destructor
    ///
    /// Destroys class.
    virtual ~GContainer(void) {}

    /// @brief Return number of objects in container
    ///
    /// @return Number of objects in container.
    virtual int size(void) const = 0;
    
    /// @brief Checks if container is empty
    ///
    /// @return True if no objects are in container, false otherwise.
    virtual bool isempty(void) const = 0;
    
    /// @brief Remove object from container
    ///
    /// @param[in] index Index.
    ///
    /// Removes the object with the specified @p index from the container.
    virtual void remove(const int& index) = 0;

    /// @brief Reserves space in the container
    ///
    /// @param[in] num Number of objects.
    ///
    /// Reserves space for @p num objects in the container.
    virtual void reserve(const int& num) = 0;
};

#endif /* GCONTAINER_HPP */
