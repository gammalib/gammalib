/***************************************************************************
 *                  GCTAEventBin.i  -  CTA event bin class                 *
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
 * @file GCTAEventBin.i
 * @brief CTA event bin class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventBin.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventBin
 *
 * @brief CTA event bin class
 ***************************************************************************/
class GCTAEventBin : public GEventBin {

    // Friend classes
    friend class GCTAEventCube;

public:
    // Constructors and destructors
    GCTAEventBin(void);
    GCTAEventBin(const GCTAEventBin& bin);
    virtual ~GCTAEventBin(void);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GCTAEventBin*      clone(void) const;
    virtual double             size(void) const;
    virtual const GCTAInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);

    // Other methods
    const double&  omega(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
};


/***********************************************************************//**
 * @brief GCTAEventBin class extension
 ***************************************************************************/
%extend GCTAEventBin {
    GCTAEventBin copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GCTAEventBin type casts
 ***************************************************************************/
%inline %{
    GCTAEventBin* cast_GCTAEventBin(GEvent* event) {
        GCTAEventBin* bin = dynamic_cast<GCTAEventBin*>(event);
        if (bin == NULL)
            throw GException::bad_type("cast_GCTAEventBin(GEvent*)",
                                       "GEvent not of type GCTAEventBin");
        return bin;
    }
%}
