/***************************************************************************
 *          GCOMEventCube.i  -  COMPTEL event bin container class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file GCOMEventCube.hpp
 * @brief COMPTEL event bin container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMEventCube.hpp"
%}


/***********************************************************************//**
 * @class GCOMEventCube
 *
 * @brief COMPTEL event bin container class
 ***************************************************************************/
class GCOMEventCube : public GEventCube {
public:
    // Constructors and destructors
    GCOMEventCube(void);
    explicit GCOMEventCube(const GFilename& filename);
    explicit GCOMEventCube(const GCOMDri& dre);
    GCOMEventCube(const GCOMEventCube& cube);
    virtual ~GCOMEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMEventCube* clone(void) const;
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
    const GCOMDri& dre(void) const;
};


/***********************************************************************//**
 * @brief GCOMEventCube class extension
 ***************************************************************************/
%extend GCOMEventCube {
    GCOMEventBin* __getitem__(int index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Event bin",
                                           index, self->size());
        }
    }
    void __setitem__(int index, const GCOMEventBin& val) {
        if (index>=0 && index < self->size()) {
            *((*self)[index]) = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Event bin",
                                           index, self->size());
        }
    }
    GCOMEventCube copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        fits = gammalib.GFits()
        self.write(fits)
        state = (fits,)
        return state
    def __setstate__(self, state):
        self.__init__()
        if not state[0].is_empty():
            self.read(state[0])
}
};
