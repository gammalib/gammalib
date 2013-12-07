/***************************************************************************
 *       GCTAModelBackground - generic CTA background model class          *
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
 * @file GCTAModelBackground.i
 * @brief CTA background model class interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelBackground.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
* @class GCTAModelBackground
*
* @brief CTA background model class
***************************************************************************/
class GCTAModelBackground : public GModelData {
    public:
    // Constructors and destructors
    GCTAModelBackground(void);
    explicit GCTAModelBackground(const GXmlElement& xml);
    explicit GCTAModelBackground(const GModelSpatial& spatial,
                                       const GModelSpectral& spectral);
    GCTAModelBackground(const GCTAModelBackground& model);
    virtual ~GCTAModelBackground(void);
    
    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GCTAModelBackground* clone(void) const;
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
    GModelSpatial* spatial(void)   const;
    GModelSpectral*  spectral(void) const;
    GModelTemporal*  temporal(void) const;
};

/***********************************************************************//**
* @brief GCTAModelBackground class extension
***************************************************************************/
%extend GCTAModelBackground {
};
