/***************************************************************************
 *           GCTAModelRadial.i - Abstract radial model base class          *
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
 * @file GCTAModelRadial.i
 * @brief Abstract radial model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadial.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadial
 *
 * @brief Abstract radial acceptance model class
 ***************************************************************************/
class GCTAModelRadial : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelRadial(void);
    GCTAModelRadial(const GCTAModelRadial& model);
    virtual ~GCTAModelRadial(void);

    // Implemented virtual methods
    virtual double           eval(const GCTAInstDir& dir,
                                  const GEnergy&     energy,
                                  const GTime&       time,
                                  const bool&        gradients = false) const;
    virtual GCTAInstDir      mc(const GEnergy&         energy,
                                const GTime&           time,
                                const GCTAObservation& roi,
                                GRan&                  ran) const;

    // Pure virtual methods
    virtual void             clear(void) = 0;
    virtual GCTAModelRadial* clone(void) const = 0;
    virtual std::string      classname(void) const = 0;
    virtual std::string      type(void) const = 0;
    virtual double           omega(void) const = 0;
    virtual double           mc_max_value(const GCTAObservation& obs) const = 0;
    virtual void             read(const GXmlElement& xml) = 0;
    virtual void             write(GXmlElement& xml) const = 0;

    // Derived class pure virtual methods
    virtual double           eval(const double& offset,
                                  const bool&   gradients = false) const = 0;
    virtual GCTAInstDir      mc(GRan& ran) const = 0;
};


/***********************************************************************//**
 * @brief GCTAModelRadial class extension
 ***************************************************************************/
%extend GCTAModelRadial {
};
