/***************************************************************************
 *                       GBase.i - GammaLib base class                     *
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
 * @file GBase.i
 * @brief Definition of interface for all GammaLib classes
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GBase.hpp"
#include "GTools.hpp"
%}


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
 * clear() Sets the object to a clean initial state
 *
 * clone() Creates a deep copy of the object
 *
 * print() Print content of object
 ***************************************************************************/
class GBase {
public:
    // Constructors and destructors
    virtual ~GBase(void);
 
    // Methods
    virtual void        clear(void) = 0;
    virtual GBase*      clone(void) const = 0;
    virtual std::string print(void) const = 0;
};


/***********************************************************************//**
 * @brief GBase class extension
 ***************************************************************************/
%extend GBase {
    char *__str__() {
        return tochar(self->print());
    }
    GBase copy() {
        return (*self);
    }
};
