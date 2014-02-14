/***************************************************************************
 *                         GEnergy.i - Energy class                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GEnergy.i
 * @brief GEnergy class python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEnergy.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GEnergy
 *
 * @brief Class that handles energies in a unit independent way.
 *
 * The GEnergy class stores an energy value in units of MeV and implements
 * methods that provide automatic conversion of the energy values in other
 * units. This makes instrument specific implementations more robust and
 * reduces the risk of unit errors.
 ***************************************************************************/
class GEnergy : public GBase {
public:
    // Constructors and destructors
    GEnergy(void);
    GEnergy(const GEnergy& eng);
    explicit GEnergy(const double& eng, const std::string& unit);
    virtual ~GEnergy(void);
 
    // Operators
    GEnergy& operator+=(const GEnergy& eng);
    GEnergy& operator-=(const GEnergy& eng);

    // Methods
    void     clear(void);
    GEnergy* clone(void) const;
    double   erg(void) const;
    double   keV(void) const;
    double   MeV(void) const;
    double   GeV(void) const;
    double   TeV(void) const;
    double   log10erg(void) const;
    double   log10keV(void) const;
    double   log10MeV(void) const;
    double   log10GeV(void) const;
    double   log10TeV(void) const;
    void     erg(const double& eng);
    void     keV(const double& eng);
    void     MeV(const double& eng);
    void     GeV(const double& eng);
    void     TeV(const double& eng);
    void     log10erg(const double& eng);
    void     log10keV(const double& eng);
    void     log10MeV(const double& eng);
    void     log10GeV(const double& eng);
    void     log10TeV(const double& eng);
    void     log10(const double& eng, const std::string& unit);
};


/***********************************************************************//**
 * @brief GEnergy class extension
 ***************************************************************************/
%extend GEnergy {
    GEnergy __add__(const GEnergy& eng) const {
        return ((*self) + eng);
    }
    GEnergy __sub__(const GEnergy& eng) const {
        return ((*self) - eng);
    }
    GEnergy __mul__(const double& factor) const {
        return ((*self) * factor);
    }
    GEnergy __div__(const double& factor) const {
        return ((*self) / factor);
    }
    bool __eq__(const GEnergy& eng) const {
        return ((*self) == eng);
    }
    bool __ne__(const GEnergy& eng) const {
        return ((*self) != eng);
    }
    bool __lt__(const GEnergy& eng) const {
        return ((*self) < eng);
    }
    bool __gt__(const GEnergy& eng) const {
        return ((*self) > eng);
    }
    bool __lte__(const GEnergy& eng) const {
        return ((*self) <= eng);
    }
    bool __gte__(const GEnergy& eng) const {
        return ((*self) >= eng);
    }
    GEnergy copy() {
        return (*self);
    }
};
