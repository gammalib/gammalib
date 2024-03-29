/***************************************************************************
 *               GCOMOad.i - COMPTEL Orbit Aspect Data class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2023 by Juergen Knodlseder                          *
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
 * @file GCOMOad.i
 * @brief COMPTEL Orbit Aspect Data class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMOad.hpp"
%}


/***********************************************************************//**
 * @class GCOMOad
 *
 * @brief COMPTEL Orbit Aspect Data class
 ***************************************************************************/
class GCOMOad : public GBase {

public:
    // Constructors and destructors
    GCOMOad(void);
    GCOMOad(const GCOMOad& oad);
    virtual ~GCOMOad(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMOad*    clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    const GTime&   tstart(void) const;
    void           tstart(const GTime& tstart);
    const GTime&   tstop(void) const;
    void           tstop(const GTime& tstop);
    const int&     tjd(void) const;
    void           tjd(const int& tjd);
    const int&     tics(void) const;
    void           tics(const int& tics);
    const float&   gcaz(void) const;
    void           gcaz(const float& gcaz);
    const float&   gcel(void) const;
    void           gcel(const float& gcel);
    const float&   georad(void) const;
    void           georad(const float& georad);
    const float&   ehora(void) const;
    void           ehora(const float& ehora);
    const GSkyDir& zaxis(void) const;
    void           zaxis(const GSkyDir& zaxis);
    const GSkyDir& xaxis(void) const;
    void           xaxis(const GSkyDir& xaxis);
    const GVector& pos(void) const;
    void           pos(const GVector& pos);
    double         theta(const GSkyDir& sky) const;
    double         phi(const GSkyDir& sky) const;
};


/***********************************************************************//**
 * @brief GCOMOad class extension
 ***************************************************************************/
%extend GCOMOad {
    GCOMOad copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.tstart(), self.tstop(), self.zaxis(), self.xaxis(),
                 self.tjd(), self.tics(), self.gcaz(), self.gcel(),
                 self.pos(), self.georad(), self.ehora())
        return state
    def __setstate__(self, state):
        self.__init__()
        self.tstart(state[0])
        self.tstop(state[1])
        self.zaxis(state[2])
        self.xaxis(state[3])
        self.tjd(state[4])
        self.tics(state[5])
        self.gcaz(state[6])
        self.gcel(state[7])
        self.pos(state[8])
        self.georad(state[9])
        self.ehora(state[10])
}
};
