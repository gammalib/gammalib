/***************************************************************************
 *              GRegistry.hpp - Interface class for registries             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GRegistry.hpp
 * @brief Interface class definition for registries
 * @author Juergen Knoedlseder
 */

#ifndef GREGISTRY_HPP
#define GREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GTypemaps.hpp"


/***********************************************************************//**
 * @class GRegistryPointer
 *
 * @brief Smart pointer for registry classes
 ***************************************************************************/
template <class T>
class GRegistryPointer {

public:
    // Constructors and destructors
    GRegistryPointer(void) : m_ptr(0) { }
    GRegistryPointer(const GRegistryPointer<T>& ptr) : m_ptr(ptr.m_ptr) { }
    virtual ~GRegistryPointer(void) { free_members(); }

    // Operators
    GRegistryPointer<T>& operator=(const GRegistryPointer<T>& ptr) {
        if (this != &ptr) {
            free_members();
            m_ptr = ptr.m_ptr;
        }
        return *this;
    }
    T& operator*() const { return *m_ptr; }
    T* operator->() const { return m_ptr; }
    T& operator[](const int& index) const { return m_ptr[index]; }

    // Methods
    void assign(T* ptr) {
        free_members();
        m_ptr = ptr;
    }

private:
    // Protected methods
    void free_members(void) {
        if (m_ptr != NULL) delete [] m_ptr;
    }

    // Private members
    T* m_ptr;  //!< Pointer
};


/***********************************************************************//**
 * @class GRegistry
 *
 * @brief Interface class for registries
 *
 * This class defines the interface for registries. A registry is a
 * container class that contains instance of derived classes. For example,
 * if three different derived classes exist for a given base class, the
 * registry will contain one instance for each derived class. Using the
 * clone mechanism, the registry may thus provide a new instance of a given
 * derived class, depending on the type of the derived class.
 *
 * The interface class requires the implementation of the following methods
 * for all registry classes:
 * 
 * size(void) Number of elements in registry class
 *
 * name(const int& index) Name of registered class by index
 *
 * print() Print content of registry
 ***************************************************************************/
class GRegistry {

public:
    /// @brief Destructor
    ///
    /// Destroys class.
    virtual ~GRegistry(void) {}

    /// @brief Return class name
    ///
    /// @return String containing the class name.
    ///
    /// Returns the class name for non-abstract classes in a human readable
    /// way.
    virtual std::string classname(void) const = 0;

    /// @brief Return number of classes in registry
    ///
    /// Returns the number of classes that have been registered in the
    /// registry.
    virtual int size(void) const = 0;

    /// @brief Return name of registered class by index
    ///
    /// @param[in] index Class index [0,...,size()-1]
    ///
    /// Returns the name of a registered class, specified by its index in
    /// the registry.
    virtual std::string name(const int& index) const = 0;

    /// @brief Print content of object
    ///
    /// @return String containing the content of the object.
    ///
    /// Formats the content in a standard way and puts this content in a
    /// C++ string that is returned.
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Implement methods
    std::string content(void) const;
};

/* __ Prototypes _________________________________________________________ */
std::ostream& operator<<(std::ostream& os, const GRegistry& registry);
GLog&         operator<<(GLog& log,        const GRegistry& registry);

#endif /* GREGISTRY_HPP */
