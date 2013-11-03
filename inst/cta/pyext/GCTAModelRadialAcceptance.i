/***************************************************************************
 *       GCTAModelRadialAcceptance.i - Radial acceptance model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialAcceptance.i
 * @brief Radial acceptance model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialAcceptance.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialAcceptance
 *
 * @brief Radial acceptance model class
 ***************************************************************************/
class GCTAModelRadialAcceptance : public GModelData {

public:
    // Constructors and destructors
    GCTAModelRadialAcceptance(void);
    explicit GCTAModelRadialAcceptance(const GXmlElement& xml);
    explicit GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                                       const GModelSpectral& spectral);
    GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model);
    virtual ~GCTAModelRadialAcceptance(void);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GCTAModelRadialAcceptance* clone(void) const;
    virtual std::string                type(void) const;
    virtual bool                       isconstant(void) const;
    virtual double                     eval(const GEvent& event,
                                            const GObservation& obs) const;
    virtual double                     eval_gradients(const GEvent& event,
                                                      const GObservation& obs) const;
    virtual double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                             const GObservation& obs) const;
    virtual GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;

    // Other methods
    GCTAModelRadial* radial(void)   const;
    GModelSpectral*  spectral(void) const;
    GModelTemporal*  temporal(void) const;
};


/***********************************************************************//**
 * @brief GCTAModelRadial class extension
 ***************************************************************************/
%extend GCTAModelRadial {
};
