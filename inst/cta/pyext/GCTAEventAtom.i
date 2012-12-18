/***************************************************************************
 *                 GCTAEventAtom.i  -  CTA event atom class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @brief CTA event bin class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventAtom.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventAtom
 *
 * @brief CTA event atom class Python interface
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
    const GCTAInstDir& dir(void) const { return m_dir; }
    const GEnergy&     energy(void) const { return m_energy; }
    const GTime&       time(void) const { return m_time; }
    void               dir(const GCTAInstDir& dir) { m_dir=dir; }
    void               energy(const GEnergy& energy) { m_energy=energy; }
    void               time(const GTime& time) { m_time=time; }
};


/***********************************************************************//**
 * @brief GCTAEventAtom class extension
 ***************************************************************************/
%extend GCTAEventAtom {
    GCTAEventAtom copy() {
        return (*self);
    }
};
