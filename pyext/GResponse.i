/***************************************************************************
 *                GResponse.i - Abstract response base class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @brief Abstract response base class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponse.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract response base class
 ***************************************************************************/
class GResponse : public GBase {

public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Pure virtual methods
    virtual void       clear(void) = 0;
    virtual GResponse* clone(void) const = 0;
    virtual bool       has_edisp(void) const = 0;
    virtual bool       has_tdisp(void) const = 0;
    virtual double     irf(const GEvent&       event,
                           const GPhoton&      photon,
                           const GObservation& obs) const = 0;
    virtual double     npred(const GPhoton&      photon,
                            const GObservation& obs) const = 0;

    // Virtual methods
    virtual double irf(const GEvent&       event,
                       const GSource&      source,
                       const GObservation& obs) const;
    virtual double irf_ptsrc(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const;
    virtual double irf_radial(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const;
    virtual double irf_elliptical(const GEvent&       event,
                                  const GSource&      source,
                                  const GObservation& obs) const;
    virtual double irf_diffuse(const GEvent&       event,
                               const GSource&      source,
                               const GObservation& obs) const;
    virtual double npred(const GSource&      source,
                         const GObservation& obs) const;
    virtual double npred_ptsrc(const GSource&      source,
                               const GObservation& obs) const;
    virtual double npred_radial(const GSource&      source,
                                const GObservation& obs) const;
    virtual double npred_elliptical(const GSource&      source,
                                    const GObservation& obs) const;
    virtual double npred_diffuse(const GSource&      source,
                                 const GObservation& obs) const;
};


/***********************************************************************//**
 * @brief GResponse class extension
 ***************************************************************************/
%extend GResponse {
};
