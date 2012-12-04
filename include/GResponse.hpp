/***************************************************************************
 *              GResponse.hpp  -  Abstract response base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
#include "GEvent.hpp"
#include "GPhoton.hpp"
//
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GModelSky.hpp"
#include "GModelPointSource.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelDiffuseSource.hpp"
#include "GIntegrand.hpp"
#include "GMatrix.hpp"

/* __ Forward declarations _______________________________________________ */
class GObservation;


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract interface for the instrument response function
 *
 * The response function provides conversion between physical parameters
 * (such as source position, flux, ...) and the measured instrumental
 * parameters (such as measured energy, photon interaction, ...).
 * For a given observation, the irf method returns the instrument response
 * for a given event and source model as function of the true photon energy
 * and the true photon arrival time.
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
    virtual GResponse& operator= (const GResponse& rsp);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GResponse*  clone(void) const = 0;
    virtual bool        hasedisp(void) const = 0;
    virtual bool        hastdisp(void) const = 0;
    virtual double      irf(const GEvent&       event,
                            const GPhoton&      photon,
                            const GObservation& obs) const = 0;
    virtual double      npred(const GPhoton&      photon,
                              const GObservation& obs) const = 0;
    virtual std::string print(void) const = 0;

    // Virtual methods
    virtual double irf(const GEvent&       event,
                       const GModelSky&    model,
                       const GEnergy&      srcEng,
                       const GTime&        srcTime,
                       const GObservation& obs) const;
    virtual double irf_ptsrc(const GEvent&            event,
                             const GModelPointSource& model,
                             const GEnergy&           srcEng,
                             const GTime&             srcTime,
                             const GObservation&      obs) const;
    virtual double irf_extended(const GEvent&               event,
                                const GModelExtendedSource& model,
                                const GEnergy&              srcEng,
                                const GTime&                srcTime,
                                const GObservation&         obs) const;
    virtual double irf_diffuse(const GEvent&              event,
                               const GModelDiffuseSource& model,
                               const GEnergy&             srcEng,
                               const GTime&               srcTime,
                               const GObservation&        obs) const;
    virtual double npred(const GModelSky&    model,
                         const GEnergy&      srcEng,
                         const GTime&        srcTime,
                         const GObservation& obs) const;
    virtual double npred_ptsrc(const GModelPointSource& model,
                               const GEnergy&           srcEng,
                               const GTime&             srcTime,
                               const GObservation&      obs) const;
    virtual double npred_extended(const GModelExtendedSource& model,
                                  const GEnergy&              srcEng,
                                  const GTime&                srcTime,
                                  const GObservation&         obs) const;
    virtual double npred_diffuse(const GModelDiffuseSource& model,
                                 const GEnergy&             srcEng,
                                 const GTime&               srcTime,
                                 const GObservation&        obs) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GResponse& rsp);
    void free_members(void);

    // Npred theta integration kernel
    class npred_kern_theta : public GIntegrand {
    public:
        npred_kern_theta(const GResponse*    rsp,
                         const GModelRadial* radial,
                         const GEnergy*      srcEng,
                         const GTime*        srcTime,
                         const GObservation* obs,
                         const GMatrix*      rot) :
                         m_rsp(rsp),
                         m_radial(radial),
                         m_srcEng(srcEng),
                         m_srcTime(srcTime),
                         m_obs(obs),
                         m_rot(rot) { return; }
        double eval(double theta);
    protected:
        const GResponse*    m_rsp;           //!< Pointer to response
        const GModelRadial* m_radial;        //!< Pointer to radial spatial model
        const GEnergy*      m_srcEng;        //!< Pointer to true photon energy
        const GTime*        m_srcTime;       //!< Pointer to true photon arrival time
        const GObservation* m_obs;           //!< Pointer to observation
        const GMatrix*      m_rot;           //!< Rotation matrix
    };

    // Npred phi integration kernel
    class npred_kern_phi : public GIntegrand {
    public:
        npred_kern_phi(const GResponse*    rsp,
                       const GEnergy*      srcEng,
                       const GTime*        srcTime,
                       const GObservation* obs,
                       const GMatrix*      rot,
                       double              theta,
                       double              sin_theta) :
                       m_rsp(rsp),
                       m_srcEng(srcEng),
                       m_srcTime(srcTime),
                       m_obs(obs),
                       m_rot(rot),
                       m_theta(theta),
                       m_cos_theta(std::cos(theta)),
                       m_sin_theta(sin_theta) { return; }
        double eval(double phi);
    protected:
        const GResponse*    m_rsp;           //!< Pointer to response
        const GEnergy*      m_srcEng;        //!< Pointer to true photon energy
        const GTime*        m_srcTime;       //!< Pointer to true photon arrival time
        const GObservation* m_obs;           //!< Pointer to observation
        const GMatrix*      m_rot;           //!< Rotation matrix
        double              m_theta;         //!< Offset angle (radians)
        double              m_cos_theta;     //!< cosine of offset angle
        double              m_sin_theta;     //!< Sine of offset angle
    };

};

#endif /* GRESPONSE_HPP */
