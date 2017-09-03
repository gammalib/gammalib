/***************************************************************************
 *                GCOMEventAtom.i - COMPTEL event atom class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMEventAtom.i
 * @brief COMPTEL event atom class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMEventAtom.hpp"
%}


/***********************************************************************//**
 * @class GCOMEventAtom
 *
 * @brief COMPTEL event atom class
 ***************************************************************************/
class GCOMEventAtom : public GEventAtom {

public:
    // Constructors and destructors
    GCOMEventAtom(void);
    GCOMEventAtom(const GCOMEventAtom& atom);
    virtual ~GCOMEventAtom(void);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCOMEventAtom*     clone(void) const;
    std::string        classname(void) const;
    const GCOMInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
};


/***********************************************************************//**
 * @brief GCOMEventAtom class extension
 ***************************************************************************/
%extend GCOMEventAtom {
    GCOMEventAtom copy() {
        return (*self);
    }
};
