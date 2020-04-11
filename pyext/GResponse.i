/***************************************************************************
 *                GResponse.i - Abstract response base class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
 * @file GResponse.i
 * @brief Abstract response base class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponse.hpp"
%}


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract instrument response base class
 ***************************************************************************/
class GResponse : public GBase {

public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GResponse*  clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual bool        use_edisp(void) const = 0;
    virtual bool        use_tdisp(void) const = 0;
    virtual double      irf(const GEvent&       event,
                            const GPhoton&      photon,
                            const GObservation& obs) const = 0;
    virtual double      nroi(const GModelSky&    model,
                             const GEnergy&      obsEng,
                             const GTime&        obsTime,
                             const GObservation& obs) const = 0;
    virtual GEbounds    ebounds(const GEnergy& obsEnergy) const = 0;

    // Virtual methods
    virtual double      convolve(const GModelSky&    model,
                                 const GEvent&       event,
                                 const GObservation& obs,
                                 const bool&         grad = true) const;
    virtual double      irf_spatial(const GEvent&       event,
                                    const GSource&      source,
                                    const GObservation& obs) const;
    virtual void        remove_response_cache(const std::string& name);
};


/***********************************************************************//**
 * @brief GResponse class extension
 ***************************************************************************/
%extend GResponse {
};
