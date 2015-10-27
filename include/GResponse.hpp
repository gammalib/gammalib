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
#include "GFunction.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GPhoton;
class GSource;
class GEnergy;
class GEbounds;
class GTime;
class GObservation;
class GModelSky;


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract instrument response base class
 *
 * The response function provides conversion between physical parameters
 * (such as source position, flux, ...) and the measured instrumental
 * parameters (such as measured energy, photon interaction, ...).
 *
 * For a given observation, the irf method returns the instrument response
 * for a given event and photon. An alternative method exists that returns
 * the response for a specific source.
 *
 * The nroi method returns the spatial integral of the instrument response
 * function times the sky model over the region of interest. This method is
 * only required for unbinned analysis.
 *
 * The ebounds method returns the true energy boundaries for a specified
 * measured event energy. This method is used for computing the energy
 * dispersion.
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
    virtual double      irf(const GEvent&       event,
                            const GSource&      source,
                            const GObservation& obs) const = 0;
    virtual double      nroi(const GModelSky&    model,
                             const GEnergy&      obsEng,
                             const GTime&        obsTime,
                             const GObservation& obs) const = 0;
    virtual GEbounds    ebounds(const GEnergy& obsEng) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Virtual methods
    virtual double      convolve(const GModelSky&    model,
                                 const GEvent&       event,
                                 const GObservation& obs,
                                 const bool&         grad = true) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GResponse& rsp);
    void   free_members(void);
    double eval_prob(const GModelSky&    model,
                     const GEvent&       event,
                     const GEnergy&      srcEng,
                     const GTime&        srcTime,
                     const GObservation& obs,
                     const bool&         grad) const;

    // Protected classes
    class edisp_kern : public GFunction {
    public:
        edisp_kern(const GResponse*    parent,
                   const GObservation* obs,
                   const GModelSky*    model,
                   const GEvent*       event,
                   const GTime&        srcTime,
                   const bool&         grad) :
                   m_parent(parent),
                   m_obs(obs),
                   m_model(model),
                   m_event(event),
                   m_srcTime(srcTime),
                   m_grad(grad) { }
        double eval(const double& x);
    protected:
        const GResponse*    m_parent;  //!< Response
        const GObservation* m_obs;     //!< Observation
        const GModelSky*    m_model;   //!< Sky model
        const GEvent*       m_event;   //!< Event
        GTime               m_srcTime; //!< True arrival time
        bool                m_grad;    //!< Gradient flag
    };
};

#endif /* GRESPONSE_HPP */
