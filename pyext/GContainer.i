/***************************************************************************
 *                 GContainer.i - Container interface class                *
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
 * @file GContainer.i
 * @brief Definition of interface for container classes
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GContainer.hpp"
%}


/***********************************************************************//**
 * @class GContainer
 *
 * @brief Interface class for container classes
 ***************************************************************************/
class GContainer : public GBase {
public:
    // Constructors and destructors
    virtual ~GContainer(void) {}

    // Pure virtual methods
    virtual int  size(void) const = 0;
    virtual bool is_empty(void) const = 0;
    virtual void remove(const int& index) = 0;
    virtual void reserve(const int& num) = 0;
};


/***********************************************************************//**
 * @brief GContainer class extension
 ***************************************************************************/
%extend GContainer {
    int __len__() {
        return (self->size());
    }
};
