/***************************************************************************
 *                  GCTAEventAtom.i - CTA event atom class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAEventAtom.i
 * @brief CTA event atom class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventAtom.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventAtom
 *
 * @brief CTA event atom class
 ***************************************************************************/
class GCTAEventAtom : public GEventAtom {
public:
    // Constructors and destructors
    GCTAEventAtom(void);
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom(void);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCTAEventAtom*     clone(void) const;
    const GCTAInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    void               dir(const GCTAInstDir& dir);
    void               energy(const GEnergy& energy);
    void               time(const GTime& time);

    // Other methods
    const int&           index(void) const;
    const unsigned long& event_id(void) const;
    const unsigned long& obs_id(void)   const;
    void                 event_id(const unsigned long& id);
    void                 obs_id(const unsigned long& id);
};


/***********************************************************************//**
 * @brief GCTAEventAtom class extension
 ***************************************************************************/
%extend GCTAEventAtom {
    GCTAEventAtom copy() {
        return (*self);
    }
};
