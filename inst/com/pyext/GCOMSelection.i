/***************************************************************************
 *              GCOMSelection.i - COMPTEL selection set class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017-2022 by Juergen Knoedlseder                         *
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
 * @file GCOMSelection.i
 * @brief COMPTEL selection set class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMSelection.hpp"
%}


/***********************************************************************//**
 * @class GCOMSelection
 *
 * @brief COMPTEL selection set class
 ***************************************************************************/
class GCOMSelection : public GBase {

public:
    // Constructors and destructors
    GCOMSelection(void);
    GCOMSelection(const GCOMSelection& select);
    virtual ~GCOMSelection(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCOMSelection* clone(void) const;
    virtual std::string    classname(void) const;

    // Other methods
    void           init_statistics(void) const;
    bool           use_event(const GCOMEventAtom& event) const;
    const double&  e1_min(void) const;
    void           e1_min(const double& e1_min);
    const double&  e1_max(void) const;
    void           e1_max(const double& e1_max);
    const double&  e2_min(void) const;
    void           e2_min(const double& e2_min);
    const double&  e2_max(void) const;
    void           e2_max(const double& e2_max);
    const int&     tof_min(void) const;
    void           tof_min(const int& tof_min);
    const int&     tof_max(void) const;
    void           tof_max(const int& tof_max);
    const int&     psd_min(void) const;
    void           psd_min(const int& psd_min);
    const int&     psd_max(void) const;
    void           psd_max(const int& psd_max);
    const int&     reflag_min(void) const;
    void           reflag_min(const int& reflag_min);
    const int&     reflag_max(void) const;
    void           reflag_max(const int& reflag_max);
    const int&     vetoflag_min(void) const;
    void           vetoflag_min(const int& vetoflag_min);
    const int&     vetoflag_max(void) const;
    void           vetoflag_max(const int& vetoflag_max);
    const int&     fpmtflag(void) const;
    void           fpmtflag(const int& fpmtflag);
    const bool&    use_d1(const int& id1) const;
    void           use_d1(const int& id1, const bool& use);
    const bool&    use_d2(const int& id2) const;
    void           use_d2(const int& id2, const bool& use);
    const GPhases& orbital_phases(void) const;
    void           orbital_phases(const GPhases& phases);
    void           orbital_period(const double& period, const GTime& time);
    double         orbital_phase(const GTime& time) const;
    const GPhases& pulsar_phases(void) const;
    void           pulsar_phases(const GPhases& phases);
    void           read(const GFitsHDU& hdu);
    void           write(GFitsHDU& hdu) const;
};


/***********************************************************************//**
 * @brief GCOMSelection class extension
 ***************************************************************************/
%extend GCOMSelection {
    GCOMSelection copy() {
        return (*self);
    }
};
