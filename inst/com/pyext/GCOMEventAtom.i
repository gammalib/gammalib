/***************************************************************************
 *                GCOMEventAtom.i - COMPTEL event atom class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
    void               clear(void);
    GCOMEventAtom*     clone(void) const;
    std::string        classname(void) const;
    const GCOMInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;

    // Other methods
    void         dir(const GCOMInstDir& dir);
    void         energy(const GEnergy& energy);
    void         time(const GTime& time);
    void         time(const int& tjd, const int& tics);
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
};
