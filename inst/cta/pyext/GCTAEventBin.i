/***************************************************************************
 *                   GCTAEventBin.i - CTA event bin class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAEventBin.i
 * @brief CTA event bin class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventBin.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventBin
 *
 * @brief CTA event bin class
 ***************************************************************************/
class GCTAEventBin : public GEventBin {

    // Friend classes
    friend class GCTAEventCube;

public:
    // Constructors and destructors
    GCTAEventBin(void);
    GCTAEventBin(const GCTAEventBin& bin);
    virtual ~GCTAEventBin(void);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GCTAEventBin*      clone(void) const;
    virtual std::string        classname(void) const;
    virtual double             size(void) const;
    virtual const GCTAInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);

    // Other methods
    const int&     ipix(void) const;
    const int&     ieng(void) const;
    const double&  solidangle(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
    const double&  weight(void) const;
    void           dir(const GCTAInstDir& dir);
    void           energy(const GEnergy& energy);
    void           time(const GTime& time);
    void           ipix(const int& ipix);
    void           ieng(const int& ieng);
    void           solidangle(const double& solidangle);
    void           ewidth(const GEnergy& ewidth);
    void           ontime(const double& ontime);
    void           weight(const double& weight);
};


/***********************************************************************//**
 * @brief GCTAEventBin class extension
 ***************************************************************************/
%extend GCTAEventBin {
    GCTAEventBin copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.ipix(), self.dir(), self.energy(), self.time(),
                 self.counts(), self.solidangle(), self.ewidth(),
                 self.ontime(), self.weight())
        return state
    def __setstate__(self, state):
        self.__init__()
        if state[0] == -1:
            self.dir(state[1])
            self.energy(state[2])
            self.time(state[3])
            self.counts(state[4])
            self.solidangle(state[5])
            self.ewidth(state[6])
            self.ontime(state[7])
            self.weight(state[8])
}
};
