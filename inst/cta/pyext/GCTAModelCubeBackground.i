/***************************************************************************
 *        GCTAModelCubeBackground - CTA cube background model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAModelCubeBackground.i
 * @brief CTA cube background model class interface definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelCubeBackground.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
* @class GCTAModelCubeBackground
*
* @brief CTA cube background model class
***************************************************************************/
class GCTAModelCubeBackground : public GModelData {
public:
    // Constructors and destructors
    GCTAModelCubeBackground(void);
    explicit GCTAModelCubeBackground(const GXmlElement& xml);
    GCTAModelCubeBackground(const GModelSpatial& spatial,
                            const GModelSpectral& spectral);
    GCTAModelCubeBackground(const GCTAModelCubeBackground& model);
    GCTAModelCubeBackground(const GCTAObservation& obs,
                            const std::string&     filename,
                            const GModelSpectral&  spectral,
                            const int&             nx_sky   = 0,
                            const int&             ny_sky   = 0,
                            const int&             n_energy = 0);
    virtual ~GCTAModelCubeBackground(void);
    
    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCTAModelCubeBackground* clone(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent& event,
                                          const GObservation& obs) const;
    virtual double                   eval_gradients(const GEvent& event,
                                                    const GObservation& obs) const;
    virtual double                   npred(const GEnergy& obsEng, const GTime& obsTime,
                                           const GObservation& obs) const;
    virtual GCTAEventList*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    
    // Other methods
    GModelSpatial*  spatial(void) const;
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;
};

/***********************************************************************//**
* @brief GCTAModelCubeBackground class extension
***************************************************************************/
%extend GCTAModelCubeBackground {
};
