/***************************************************************************
 *         GXXXEventCube.i - [INSTRUMENT] event bin container class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXEventCube.i
 * @brief [INSTRUMENT] event bin container class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXEventCube.hpp"
%}


/***********************************************************************//**
 * @class GXXXEventCube
 *
 * @brief [INSTRUMENT] event bin container class
 ***************************************************************************/
class GXXXEventCube : public GEventCube {

public:
    // Constructors and destructors
    GXXXEventCube(void);
    explicit GXXXEventCube(const GFilename& filename);
    GXXXEventCube(const GXXXEventCube& cube);
    virtual ~GXXXEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GXXXEventCube* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;

    // Other methods
    // TODO: Copy methods from GXXXEventCube.hpp file
};


/***********************************************************************//**
 * @brief GXXXEventCube class extension
 ***************************************************************************/
%extend GXXXEventCube {
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
    GXXXEventCube copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = ()
        return state
    def __setstate__(self, state):
        self.__init__()
}
};
