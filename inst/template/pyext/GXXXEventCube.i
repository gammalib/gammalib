/***************************************************************************
 *            GXXXEventCube.i  -  XXX event bin container class            *
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
 * @file GXXXEventCube.hpp
 * @brief XXX event bin container class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GException.hpp"
#include "GXXXEventCube.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventCube
 *
 * @brief CTA event bin container class Python interface
 ***************************************************************************/
class GXXXEventCube : public GEventCube {
public:
    // Constructors and destructors
    GXXXEventCube(void);
    explicit GXXXEventCube(const std::string& filename);
    GXXXEventCube(const GXXXEventCube& cube);
    virtual ~GXXXEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GXXXEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(int axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename, bool clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
};


/***********************************************************************//**
 * @brief GXXXEventCube class extension
 ***************************************************************************/
%extend GXXXEventCube {
    GXXXEventCube copy() {
        return (*self);
    }
    GXXXEventBin* __getitem__(int index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(int index, const GXXXEventBin& val) {
        if (index>=0 && index < self->size()) {
            *((*self)[index]) = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
};
