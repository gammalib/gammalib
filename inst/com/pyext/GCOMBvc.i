/***************************************************************************
 *          GCOMBvc.i - COMPTEL Solar System Barycentre Data class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvc.i
 * @brief COMPTEL Solar System Barycentre Data class definition
 * @author Juergen Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMBvc.hpp"
%}


/***********************************************************************//**
 * @class GCOMBvc
 *
 * @brief COMPTEL Solar System Barycentre Data class
 ***************************************************************************/
class GCOMBvc : public GBase {

public:
    // Constructors and destructors
    GCOMBvc(void);
    GCOMBvc(const GCOMBvc& bvc);
    virtual ~GCOMBvc(void);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMBvc*    clone(void) const;
    virtual std::string classname(void) const;

    // Other methods
    const GTime&   time(void) const;
    void           time(const GTime& time);
    const int&     tjd(void) const;
    void           tjd(const int& tjd);
    const int&     tics(void) const;
    void           tics(const int& tics);
    const GVector& ssb(void) const;
    void           ssb(const GVector& ssb);
    const double&  tdelta(void) const;
    void           tdelta(const double& tdelta);
    double         tdelta(const GSkyDir& dir) const;
};


/***********************************************************************//**
 * @brief GCOMBvc class extension
 ***************************************************************************/
%extend GCOMBvc {
    GCOMBvc copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = {'time':   self.time(),
                 'tjd':    self.tjd(),
                 'tics':   self.tics(),
                 'ssb':    self.ssb(),
                 'tdelta': self.tdelta()}
        return state
    def __setstate__(self, state):
        self.__init__()
        self.time(state['time'])
        self.tjd(state['tjd'])
        self.tics(state['tics'])
        self.ssb(state['ssb'])
        self.tdelta(state['tdelta'])
}
};
