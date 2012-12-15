/***************************************************************************
 *                    GPhotons.i - Photon container class                  *
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
 * @file GPhoton.i
 * @brief Photon container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPhotons.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPhotons
 *
 * @brief Container class for photons.
 ***************************************************************************/
class GPhotons : public GBase {
public:
    // Constructors and destructors
    GPhotons(void);
    GPhotons(const GPhotons& photons);
    virtual ~GPhotons(void);
 
    // Operators
    GPhotons&      operator=(const GPhotons& photons);
    GPhoton&       operator[](const int& index);
    const GPhoton& operator[](const int& index) const;

    // Methods
    void           clear(void);
    GPhotons*      clone(void) const;
    int            size(void) const;
    void           append(const GPhoton& photon);
    void           reserve(const int& number);
};


/***********************************************************************//**
 * @brief GPhotons class extension
 ***************************************************************************/
%extend GPhotons {
    char *__str__() {
        return tochar(self->print());
    }
    GPhoton& __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GPhoton& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GPhotons copy() {
        return (*self);
    }
};
