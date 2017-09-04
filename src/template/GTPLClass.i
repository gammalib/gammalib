/***************************************************************************
 *                        GTPLClass.i - [WHAT] class                       *
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
 * @file GTPLClass.i
 * @brief [WHAT] class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTPLClass.hpp"
%}


/***********************************************************************//**
 * @class GTPLClass
 *
 * @brief [WHAT] class
 ***************************************************************************/
class GTPLClass : public GRoi {

public:
    // Constructors and destructors
    GTPLClass(void);
    GTPLClass(const GTPLClass& TPL_OBJECT);
    virtual ~GTPLClass(void);

    // Implemented pure virtual base class methods
    virtual void       clear(void);
    virtual GTPLClass* clone(void) const;

    // Other methods
    // TODO: Copy methods from GTPLClass.hpp file
};


/***********************************************************************//**
 * @brief GTPLClass class extension
 ***************************************************************************/
%extend GTPLClass {
    GTPLClass copy() {
        return (*self);
    }
};
