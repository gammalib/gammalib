/***************************************************************************
 *               GRegistry.i - Interface class for registries              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GRegistry.i
 * @brief Interface class definition for registries
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GRegistry.hpp"
#include "GTools.hpp"
%}


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
    // Constructors and destructors
    virtual ~GRegistry(void);
 
    // Methods
    virtual int         size(void) const = 0;
    virtual std::string name(const int& index) const = 0;
};


/***********************************************************************//**
 * @brief GRegistry class extension
 ***************************************************************************/
%extend GRegistry {
    char *__str__() {
        return gammalib::tochar(self->print());
    }
    int __len__() {
        return (self->size());
    }
};
