/***************************************************************************
 *             GCTAResponse.i - CTA response abstract base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAResponse.i
 * @brief CTA response abstract base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAResponse.hpp"
%}


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief CTA instrument response function class
 ***************************************************************************/
class GCTAResponse : public GResponse {
public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    virtual ~GCTAResponse(void);

    // Pure virtual methods
    virtual void          clear(void) = 0;
    virtual GCTAResponse* clone(void) const = 0;
    virtual bool          use_edisp(void) const = 0;
    virtual bool          use_tdisp(void) const = 0;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const = 0;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const = 0;
    virtual void          read(const GXmlElement& xml) = 0;
    virtual void          write(GXmlElement& xml) const = 0;
};


/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponse {
};
