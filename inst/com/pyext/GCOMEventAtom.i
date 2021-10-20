/***************************************************************************
 *                GCOMEventAtom.i - COMPTEL event atom class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2021 by Juergen Knoedlseder                         *
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
    virtual void                 clear(void);
    virtual GCOMEventAtom*       clone(void) const;
    virtual std::string          classname(void) const;
    virtual const GCOMInstDir&   dir(void) const;
    virtual const GEnergy&       energy(void) const;
    virtual const GTime&         time(void) const;
    virtual const GPolarization& polarization(void) const;

    // Other methods
    void         dir(const GCOMInstDir& dir);
    void         energy(const GEnergy& energy);
    void         time(const GTime& time);
    void         time(const int& tjd, const int& tics);
    void         polarization(const GPolarization& polarizations);
    void         phibar(const float& phibar);
    const float& phibar(void) const;
    void         phi(const float& phi);
    const float& phi(void) const;
    void         theta(const float& theta);
    const float& theta(void) const;
    void         eha(const float& eha);
    const float& eha(void) const;
    void         e1(const float& e1);
    const float& e1(void) const;
    void         e2(const float& e2);
    const float& e2(void) const;
    void         psd(const float& psd);
    const float& psd(void) const;
    void         tof(const float& tof);
    const float& tof(void) const;
    void         modcom(const int& modcom);
    const int&   modcom(void) const;
    void         reflag(const int& reflag);
    const int&   reflag(void) const;
    void         veto(const int& veto);
    const int&   veto(void) const;
};


/***********************************************************************//**
 * @brief GCOMEventAtom class extension
 ***************************************************************************/
%extend GCOMEventAtom {
    GCOMEventAtom copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.dir(), self.energy(), self.time(), self.polarization(),
                 self.e1(), self.e2(), self.phibar(), self.theta(),
                 self.phi(), self.eha(), self.psd(), self.tof(),
                 self.modcom(), self.reflag(), self.veto())
        return state
    def __setstate__(self, state):
        self.__init__()
        self.dir(state[0])
        self.energy(state[1])
        self.time(state[2])
        self.polarization(state[3])
        self.e1(state[4])
        self.e2(state[5])
        self.phibar(state[6])
        self.theta(state[7])
        self.phi(state[8])
        self.eha(state[9])
        self.psd(state[10])
        self.tof(state[11])
        self.modcom(state[12])
        self.reflag(state[13])
        self.veto(state[14])
}
};
