/***************************************************************************
 *             GXXXRoi.i - [INSTRUMENT] region of interest class           *
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
 * @file GXXXRoi.i
 * @brief [INSTRUMENT] region of interest class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXRoi.hpp"
%}


/***********************************************************************//**
 * @class GXXXRoi
 *
 * @brief [INSTRUMENT] region of interest class
 ***************************************************************************/
class GXXXRoi : public GRoi {

public:
    // Constructors and destructors
    GXXXRoi(void);
    GXXXRoi(const GXXXRoi& roi);
    virtual ~GXXXRoi(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GXXXRoi*    clone(void) const;
    virtual std::string classname(void) const;
    virtual bool        contains(const GEvent& event) const;

    // Other methods
    // TODO: Copy methods from GXXXRoi.hpp file
};


/***********************************************************************//**
 * @brief GXXXRoi class extension
 ***************************************************************************/
%extend GXXXRoi {
    GXXXRoi copy() {
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
