/***************************************************************************
 *               GResponse.hpp - Abstract response base class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2015 by Juergen Knoedlseder                         *
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
 * @file GResponse.hpp
 * @brief Abstract response base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GRESPONSE_HPP
#define GRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GPhoton;
class GEnergy;
class GTime;
class GObservation;
class GModelSky;

// To be removed later ...
class GSource;
class GEbounds;

/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract instrument response base class
 *
 * The response function provides conversion between physical parameters
 * (such as source position, flux, ...) and the measured instrumental
 * parameters (such as measured energy, photon interaction, ...).
 * For a given observation, the irf method returns the instrument response
 * for a given event and .
 * The npred method returns the integral of the instrument response function
 * over the dataspace. This method is only required for unbinned analysis.
 ***************************************************************************/
class GResponse : public GBase {

public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Operators
    virtual GResponse& operator=(const GResponse& rsp);

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
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Virtual methods
    virtual double      convolve(const GModelSky&    model,
                                 const GEvent&       event,
                                 const GObservation& obs,
                                 const bool&         grad = true) const;

    // Old methods that will become obsolete
    virtual double   irf(const GEvent&       event,
                         const GSource&      source,
                         const GObservation& obs) const;
    virtual double   irf_ptsrc(const GEvent&       event,
                               const GSource&      source,
                               const GObservation& obs) const;
    virtual double   irf_radial(const GEvent&       event,
                                const GSource&      source,
                                const GObservation& obs) const;
    virtual double   irf_elliptical(const GEvent&       event,
                                    const GSource&      source,
                                    const GObservation& obs) const;
    virtual double   irf_diffuse(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const;
    virtual GEbounds ebounds_src(const GEnergy& obsEng) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GResponse& rsp);
    void free_members(void);
};

#endif /* GRESPONSE_HPP */
