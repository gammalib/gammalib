/***************************************************************************
 *         GLATEventList.i  -  Fermi/LAT event atom container class        *
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
 * @file GLATEventList.i
 * @brief Fermi/LAT event atom container class interface definition
 * @author J. Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventList.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GLATEventList
 *
 * @brief Fermi/LAT event atom container class
 ***************************************************************************/
class GLATEventList : public GEventList {
public:
    // Constructors and destructors
    GLATEventList(void);
    GLATEventList(const GLATEventList& list);
    virtual ~GLATEventList(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GLATEventList* clone(void) const;
    virtual int            size(void) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename, bool clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual void           roi(const GRoi& roi);
    virtual const GLATRoi& roi(void) const;
};


/***********************************************************************//**
 * @brief GLATEventList class extension
 ***************************************************************************/
%extend GLATEventList {
    char *__str__() {
        return tochar(self->print());
    }
    GLATEventList copy() {
        return (*self);
    }
    GLATEventAtom* __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return ((*self)[index]);
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const GLATEventAtom& value) {
        if (index>=0 && index < self->size()) {
            *((*self)[index]) = value;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
};


/***********************************************************************//**
 * @brief GLATEventList type casts
 ***************************************************************************/
%inline %{
    GLATEventList* cast_GLATEventList(GEvents* events) {
        GLATEventList* list = dynamic_cast<GLATEventList*>(events);
        if (list == NULL) {
            throw GException::bad_type("cast_GLATEventList(GEvents*)",
                                       "GEvents not of type GLATEventList");
        }
        return list;
    }
%}
