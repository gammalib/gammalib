/***************************************************************************
 *                      GBase.hpp - GammaLib base class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GBase.hpp
 * @brief Definition of interface for all GammaLib classes
 * @author Juergen Knoedlseder
 */

#ifndef GBASE_HPP
#define GBASE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GTypemaps.hpp"


/***********************************************************************//**
 * @class GBase
 *
 * @brief Interface class for all GammaLib classes
 *
 * This class defines the interface for all GammaLib classes. It is an
 * abstract base class from which all other GammaLib classes will be
 * derived. The interface class imposes on all GammaLib classes to
 * implement the following methods:
 *
 *     clear      - Clear container
 *     clone      - Clones container
 *     class_name - Returns the class name
 *     print      - Print container content
 *
 ***************************************************************************/
class GBase {

public:
    /// @brief Destructor
    ///
    /// Destroys class.
    virtual ~GBase(void) {}
 
    /// @brief Clear object
    ///
    /// Sets the object to a clean initial state. After calling the method
    /// the object will be in the same state as it were if an empty instance
    /// of the object would have been created.
    virtual void        clear(void) = 0;

    /// @brief Clones object
    ///
    /// @return Pointer to deep copy of object.
    ///
    /// Creates a deep copy of the object and returns a pointer to the
    /// object.
    virtual GBase*      clone(void) const = 0;

    /// @brief Return class name
    ///
    /// @return String containing the class name.
    ///
    /// Returns the class name for non-abstract classes in a human readable
    /// way.
    virtual std::string classname(void) const = 0;

    /// @brief Print content of object
    ///
    /// @param[in] chatter Chattiness (defaults to NORMAL).
    /// @return String containing the content of the object.
    ///
    /// Formats the content in a standard way and puts this content in a
    /// C++ string that is returned.
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;
};

/* __ Forward declarations _______________________________________________ */
class GLog;

/* __ Prototypes _________________________________________________________ */
std::ostream& operator<<(std::ostream& os, const GBase& base);
GLog&         operator<<(GLog& log,        const GBase& base);

#endif /* GBASE_HPP */
