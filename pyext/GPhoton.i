/***************************************************************************
 *                   GPhoton.i  -  Photon class python I/F                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
 * @brief GPhoton class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPhoton.hpp"
#include "GTools.hpp"
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
class GPhoton {
public:
    // Constructors and destructors
    GPhoton(void);
    GPhoton(const GPhoton& ph);
    virtual ~GPhoton(void);
 
    // Methods
    void           clear(void);
    const GSkyDir& dir(void) const { return m_dir; }
    const GEnergy& energy(void) const { return m_energy; }
    const GTime&   time(void) const { return m_time; }
    void           dir(const GSkyDir& dir) { m_dir=dir; }
    void           energy(const GEnergy& energy) { m_energy=energy; }
    void           time(const GTime& time) { m_time=time; }
};


/***********************************************************************//**
 * @brief GPhoton class extension
 ***************************************************************************/
%extend GPhoton {
    char *__str__() {
        return tochar(self->print());
    }
    bool __eq__(const GPhoton& ph) const {
        return ((*self) == ph);
    }
    bool __ne__(const GPhoton& ph) const {
        return ((*self) != ph);
    }
    GPhoton copy() {
        return (*self);
    }
};


/***************************************************************************
 *                                 Typedefs                                *
 ***************************************************************************/
typedef std::vector<GPhoton> GPhotons;
%template(GPhotons) std::vector<GPhoton>;
