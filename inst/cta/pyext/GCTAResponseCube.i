/***************************************************************************
 *       GCTAResponseCube.i - CTA cube-style response function class       *
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
 * @file GCTAResponseCube.i
 * @brief CTA cube-style response function class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAResponseCube.hpp"
%}


/***********************************************************************//**
 * @class GCTAResponseCube
 *
 * @brief CTA cube-style response function class
 ***************************************************************************/
class GCTAResponseCube : public GCTAResponse {
public:
    // Constructors and destructors
    GCTAResponseCube(void);
    GCTAResponseCube(const GCTAResponseCube& rsp);
    virtual ~GCTAResponseCube(void);

    // Implement pure virtual base class methods
    virtual void              clear(void);
    virtual GCTAResponseCube* clone(void) const;
    virtual bool              use_edisp(void) const;
    virtual bool              use_tdisp(void) const;
    virtual double            irf(const GEvent&       event,
                                  const GPhoton&      photon,
                                  const GObservation& obs) const;
    virtual double            npred(const GPhoton&      photon,
                                    const GObservation& obs) const;
    virtual void              read(const GXmlElement& xml);
    virtual void              write(GXmlElement& xml) const;

    // Other Methods
    const GCTAExposure&       exposure(void) const;
    void                      exposure(const GCTAExposure& exposure);
    const GCTAMeanPsf&        psf(void) const;
    void                      psf(const GCTAMeanPsf& psf);
};


/***********************************************************************//**
 * @brief GCTAResponseCube class extension
 ***************************************************************************/
%extend GCTAResponseCube {
    GCTAResponseCube copy() {
        return (*self);
    }
};
