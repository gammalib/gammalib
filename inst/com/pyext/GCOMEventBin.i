/***************************************************************************
 *                GCOMEventBin.i  -  COMPTEL event bin class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
 * @file GCOMEventBin.hpp
 * @brief COMPTEL event bin class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMEventBin.hpp"
%}


/***********************************************************************//**
 * @class GCOMEventBin
 *
 * @brief COMPTEL event bin class
 ***************************************************************************/
class GCOMEventBin : public GEventBin {
public:
    // Constructors and destructors
    GCOMEventBin(void);
    GCOMEventBin(const GCOMEventBin& bin);
    virtual ~GCOMEventBin(void);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GCOMEventBin*      clone(void) const;
    virtual std::string        classname(void) const;
    virtual double             size(void) const;
    virtual const GCOMInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);

    // Other methods
    const int&     index(void) const;
    const double&  solidangle(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
    void           index(const int& index);
    void           dir(const GCOMInstDir& dir);
    void           energy(const GEnergy& energy);
    void           time(const GTime& time);
    void           solidangle(const double& solidangle);
    void           ewidth(const GEnergy& ewidth);
    void           ontime(const double& ontime);
};


/***********************************************************************//**
 * @brief GCOMEventBin class extension
 ***************************************************************************/
%extend GCOMEventBin {
    GCOMEventBin copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.index(), self.dir(), self.energy(), self.time(),
                 self.counts(), self.solidangle(), self.ewidth(),
                 self.ontime())
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
}
};
