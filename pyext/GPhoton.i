/***************************************************************************
 *                         GPhoton.i - Photon class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GPhoton.i
 * @brief Photon class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPhoton.hpp"
%}


/***********************************************************************//**
 * @class GPhoton
 *
 * @brief Class that handles photons.
 *
 * The GPhoton class stores the physical attributes of a photon such as the
 * photon arrival direction, its energy and its arrival time. This class is
 * mainly used for Monte Carlo simulations.
 ***************************************************************************/
class GPhoton : public GBase {
public:
    // Constructors and destructors
    GPhoton(void);
    GPhoton(const GSkyDir& dir, const GEnergy& energy, const GTime& time,
            const int& mc_id = -1);
    GPhoton(const GPhoton& photon);
    virtual ~GPhoton(void);
 
    // Methods
    void           clear(void);
    GPhoton*       clone(void) const;
    std::string    classname(void) const;
    const GSkyDir& dir(void) const;
    const GEnergy& energy(void) const;
    const GTime&   time(void) const;
    const int&     mc_id(void) const;
    void           dir(const GSkyDir& dir);
    void           energy(const GEnergy& energy);
    void           time(const GTime& time);
    void           mc_id(const int& mc_id);
};


/***********************************************************************//**
 * @brief GPhoton class extension
 ***************************************************************************/
%extend GPhoton {
    bool __eq__(const GPhoton& photon) const {
        return ((*self) == photon);
    }
    bool __ne__(const GPhoton& photon) const {
        return ((*self) != photon);
    }
    GPhoton copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        args = self.dir(), self.energy(), self.time(), self.mc_id()
        return args
    def __setstate__(self, state):
        self.__init__()
        self.dir(state[0])
        self.energy(state[1])
        self.time(state[2])
        self.mc_id(state[3])
}
};
