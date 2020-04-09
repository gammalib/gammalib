/***************************************************************************
 *         GXXXInstDir.i - [INSTRUMENT] instrument direction class         *
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
 * @file GXXXInstDir.i
 * @brief [INSTRUMENT] instrument direction class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXInstDir.hpp"
%}


/***********************************************************************//**
 * @class GXXXInstDir
 *
 * @brief [INSTRUMENT] instrument direction class
 ***************************************************************************/
class GXXXInstDir : public GInstDir {

public:
    // Constructors and destructors
    GXXXInstDir(void);
    GXXXInstDir(const GXXXInstDir& dir);
    virtual ~GXXXInstDir(void);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GXXXInstDir* clone(void) const;
    virtual double       hash(void) const;
    virtual std::string  classname(void) const;

    // Other methods
    // TODO: Copy methods from GXXXInstDir.hpp file
};


/***********************************************************************//**
 * @brief GXXXInstDir class extension
 ***************************************************************************/
%extend GXXXInstDir {
    GXXXInstDir copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.classname()) # TODO: Replace by appropriate class members
        return state
    def __setstate__(self, state):
        self.__init__()
}
};
