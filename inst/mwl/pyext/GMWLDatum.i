/***************************************************************************
 *          GMWLDatum.i - Multi-wavelength spectral point class            *
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
 * @file GMWLDatum.i
 * @brief Multi-wavelength spectral point class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLDatum.hpp"
%}


/***********************************************************************//**
 * @class GMWLDatum
 *
 * @brief Multi-wavelength spectral point class
 *
 * This class defines a spectral data point for the multi-wavelength
 * interface. It derives from the abstract GEventBin base class.
 ***************************************************************************/
class GMWLDatum : public GEventBin {
public:
    // Constructors and destructors
    GMWLDatum(void);
    GMWLDatum(const GMWLDatum& datum);
    virtual ~GMWLDatum(void);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GMWLDatum*         clone(void) const;
    virtual double             size(void) const;
    virtual const GInstDir&    dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);

    // Other methods
    const GEnergy& energy_err(void) const;
    const double&  flux(void) const;
    const double&  flux_err(void) const;
    void           flux(const double& flux);
    void           flux_err(const double& flux_err);
};


/***********************************************************************//**
 * @brief GMWLDatum class extension
 ***************************************************************************/
%extend GMWLDatum {
    GMWLDatum copy() {
        return (*self);
    }
};
