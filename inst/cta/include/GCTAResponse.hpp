/***************************************************************************
 *            GCTAResponse.hpp - CTA response abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2019 by Juergen Knoedlseder                         *
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
 * @file GCTAResponse.hpp
 * @brief CTA response abstract base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSE_HPP
#define GCTARESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GResponse.hpp"
#include "GCTAResponseCache.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSky;
class GEnergy;
class GTime;
class GEvent;
class GPhoton;
class GSource;
class GObservation;
class GXmlElement;


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief CTA instrument response function class
 *
 * This class defines the interface for the CTA response function. It
 * provides an abstract base class to CTA response function classes.
 ***************************************************************************/
class GCTAResponse : public GResponse {

public:
    // Constructors and destructors
    GCTAResponse(void);
    GCTAResponse(const GCTAResponse& rsp);
    virtual ~GCTAResponse(void);

    // Operators
    virtual GCTAResponse& operator=(const GCTAResponse & rsp);

    // Pure virtual methods
    virtual void          clear(void) = 0;
    virtual GCTAResponse* clone(void) const = 0;
    virtual std::string   classname(void) const = 0;
    virtual bool          is_valid(void) const = 0;
    virtual bool          use_edisp(void) const = 0;
    virtual bool          use_tdisp(void) const = 0;
    virtual bool          apply_edisp(void) const = 0;
    virtual void          apply_edisp(const bool& apply_edisp) const = 0;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const = 0;
    virtual double        nroi(const GModelSky&    model,
                               const GEnergy&      obsEng,
                               const GTime&        obsTime,
                               const GObservation& obs) const = 0;
    virtual GEbounds      ebounds(const GEnergy& obsEng) const = 0;
    virtual void          read(const GXmlElement& xml) = 0;
    virtual void          write(GXmlElement& xml) const = 0;
    virtual std::string   print(const GChatter& chatter = NORMAL) const = 0;

    // Implemeted methods
    virtual double irf(const GEvent&       event,
                       const GSource&      source,
                       const GObservation& obs) const;
    virtual void   remove_response_cache(const std::string& name);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAResponse& rsp);
    void free_members(void);

    // Pure virtual protected methods
    virtual double irf_ptsrc(const GEvent&       event,
                             const GSource&      source,
                             const GObservation& obs) const = 0;
    virtual double irf_radial(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const = 0;
    virtual double irf_elliptical(const GEvent&       event,
                                  const GSource&      source,
                                  const GObservation& obs) const = 0;
    virtual double irf_diffuse(const GEvent&       event,
                               const GSource&      source,
                               const GObservation& obs) const = 0;

    // Implemented protected methods
    virtual double irf_composite(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const;

    // Protected members
    bool m_use_irf_cache;   //!< Control usage of irf cache
    bool m_use_nroi_cache;  //!< Control usage of nroi cache

    // Cache for irf(GEvent&, GSource&, GObservation&) and
    // nroi(GModelSky&, GEnergy&, GTime&, GEnergy&, GTime&, GObservation&)
    // computations
    mutable GCTAResponseCache m_irf_cache;
    mutable GCTAResponseCache m_nroi_cache;
};

#endif /* GCTARESPONSE_HPP */
