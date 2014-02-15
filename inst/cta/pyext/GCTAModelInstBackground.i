/***************************************************************************
 *    GCTAModelInstBackground.i - CTA instrument background model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTAModelInstBackground.i
 * @brief CTA instrument background model class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelInstBackground.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelInstBackground
 *
 * @brief CTA instrument background model class
 ***************************************************************************/
class GCTAModelInstBackground : public GModelData {
public:
    // Constructors and destructors
    GCTAModelInstBackground(void);
    explicit GCTAModelInstBackground(const GXmlElement& xml);
    explicit GCTAModelInstBackground(const GModelSpectral& spectral);
    GCTAModelInstBackground(const GCTAModelInstBackground& bgd);
    virtual ~GCTAModelInstBackground(void);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCTAModelInstBackground* clone(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent& event,
                                          const GObservation& obs) const;
    virtual double                   eval_gradients(const GEvent& event,
                                                    const GObservation& obs) const;
    virtual double                   npred(const GEnergy& obsEng,
                                           const GTime& obsTime,
                                           const GObservation& obs) const;
    virtual GCTAEventList*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;

    // Other methods
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;
};

/***********************************************************************//**
* @brief GCTAModelInstBackground class extension
***************************************************************************/
%extend GCTAModelInstBackground {
    GCTAModelInstBackground copy() {
        return (*self);
    }
};
