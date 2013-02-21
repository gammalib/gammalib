/***************************************************************************
 *                      GTimes.i - Time container class                    *
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
 * @file GTimes.i
 * @brief Time container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTimes.hpp"
#include "GException.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GTimes
 *
 * @brief Container class for times.
 ***************************************************************************/
class GTimes : public GContainer {

public:
    // Constructors and destructors
    GTimes(void);
    GTimes(const GTimes& times);
    virtual ~GTimes(void);
 
    // Methods
    void    clear(void);
    GTimes* clone(void) const;
    int     size(void) const;
    bool    isempty(void) const;
    void    append(const GTime& time);
    void    insert(const int& index, const GTime& time);
    void    remove(const int& index);
    void    reserve(const int& num);
    void    extend(const GTimes& times);
};


/***********************************************************************//**
 * @brief GTimes class extension
 ***************************************************************************/
%extend GTimes {
    char *__str__() {
        return tochar(self->print());
    }
    GTime& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GTime& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    int __len__() {
        return (self->size());
    }
    GTimes copy() {
        return (*self);
    }
};
