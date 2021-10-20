/***************************************************************************
 *                  GCTAEventAtom.i - CTA event atom class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2021 by Juergen Knoedlseder                         *
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
    GCTAEventAtom(const GCTAInstDir& dir, const GEnergy& energy, const GTime& time);
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom(void);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GCTAEventAtom*       clone(void) const;
    virtual std::string          classname(void) const;
    virtual const GCTAInstDir&   dir(void) const;
    virtual const GEnergy&       energy(void) const;
    virtual const GTime&         time(void) const;
    virtual const GPolarization& polarization(void) const;

    // Other methods
    const int&           index(void) const;
    const unsigned long& event_id(void) const;
    const int&           mc_id(void) const;
    const float&         phase(void) const;
    void                 index(const int& index);
    void                 event_id(const unsigned long& id);
    void                 mc_id(const int& id);
    void                 phase(const float& phase);
    void                 dir(const GCTAInstDir& dir);
    void                 energy(const GEnergy& energy);
    void                 time(const GTime& time);
    void                 polarization(const GPolarization& polarization);
};


/***********************************************************************//**
 * @brief GCTAEventAtom class extension
 ***************************************************************************/
%extend GCTAEventAtom {
    GCTAEventAtom copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        args = self.index(), self.dir(), self.energy(), self.time(), \
               self.event_id(), self.mc_id(), self.phase()
        return args
    def __setstate__(self, state):
        self.__init__()
        self.index(state[0])
        self.dir(state[1])
        self.energy(state[2])
        self.time(state[3])
        self.event_id(state[4])
        self.mc_id(state[5])
        self.phase(state[6])
}
};
