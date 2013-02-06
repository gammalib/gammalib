/***************************************************************************
 *            GSPIEventCube.i  -  SPI event bin container class            *
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
 * @file GSPIEventCube.hpp
 * @brief SPI event bin container class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GException.hpp"
#include "GSPIEventCube.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventCube
 *
 * @brief CTA event bin container class Python interface
 ***************************************************************************/
class GSPIEventCube : public GEventCube {
public:
    // Constructors and destructors
    GSPIEventCube(void);
    explicit GSPIEventCube(const std::string& filename);
    GSPIEventCube(const GSPIEventCube& cube);
    virtual ~GSPIEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GSPIEventCube* clone(void) const;
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
 * @brief GSPIEventCube class extension
 ***************************************************************************/
%extend GSPIEventCube {
    GSPIEventCube copy() {
        return (*self);
    }
    GSPIEventBin* __getitem__(int index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(int index, const GSPIEventBin& val) {
        if (index>=0 && index < self->size()) {
            *((*self)[index]) = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
};
