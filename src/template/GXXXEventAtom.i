/***************************************************************************
 *              GXXXEventAtom.i - [INSTRUMENT] event atom class            *
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
 * @file GXXXEventAtom.i
 * @brief [INSTRUMENT] event atom class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXEventAtom.hpp"
%}


/***********************************************************************//**
 * @class GXXXEventAtom
 *
 * @brief [INSTRUMENT] event atom class
 ***************************************************************************/
class GXXXEventAtom : public GEventAtom {

public:
    // Constructors and destructors
    GXXXEventAtom(void);
    GXXXEventAtom(const GXXXEventAtom& atom);
    virtual ~GXXXEventAtom(void);

    // Implemented pure virtual base class methods
    void               clear(void);
    GXXXEventAtom*     clone(void) const;
    std::string        classname(void) const;
    const GXXXInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
};


/***********************************************************************//**
 * @brief GXXXEventAtom class extension
 ***************************************************************************/
%extend GXXXEventAtom {
    GXXXEventAtom copy() {
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
