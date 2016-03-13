/***************************************************************************
 *                       GBase.i - GammaLib base class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * The GBase class defines the interface for all GammaLib classes. It is an
 * abstract base class from which most GammaLib classes will be derived. The
 * interface class imposes on all GammaLib classes the implementation of the
 * following methods:
 * 
 * clear() Sets the object to a clean initial state
 *
 * clone() Creates a deep copy of the object
 *
 * classname() Returns the mane of the class
 ***************************************************************************/
class GBase {
public:
    // Constructors and destructors
    virtual ~GBase(void);

    // Methods
    virtual void        clear(void) = 0;
    virtual GBase*      clone(void) const = 0;
    virtual std::string classname(void) const = 0;
};


/***********************************************************************//**
 * @brief GBase class extension
 ***************************************************************************/
%extend GBase {
    char *__str__() {
        return gammalib::tochar(self->print());
    }
};
