/***************************************************************************
 *               GResponse.hpp - Abstract response base class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2021 by Juergen Knoedlseder                         *
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
#include "GTime.hpp"
#include "GResponseCache.hpp"
#include "GResponseVectorCache.hpp"
#include "GFunctions.hpp"

/* __ Forward declarations _______________________________________________ */
class GVector;
class GMatrix;
class GMatrixSparse;
class GEvent;
class GPhoton;
class GSource;
class GSkyDir;
class GEnergy;
class GEbounds;
class GObservation;
class GModelPar;
class GModelSky;
class GModelSpatialRadial;
class GModelSpatialElliptical;


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
    virtual GVector     convolve(const GModelSky&    model,
                                 const GObservation& obs,
                                 GMatrixSparse*      gradients = NULL) const;
    virtual double      irf_spatial(const GEvent&       event,
                                    const GSource&      source,
                                    const GObservation& obs) const;
    virtual GVector     irf_spatial(const GModelSky&    model,
                                    const GObservation& obs,
                                    GMatrix*            gradients = NULL) const;
    virtual void        remove_response_cache(const std::string& name);

protected:
    // Protected methods
    void     init_members(void);
    void     copy_members(const GResponse& rsp);
    void     free_members(void);
    double   eval_prob(const GModelSky&    model,
                       const GEvent&       event,
                       const GEnergy&      srcEng,
                       const GTime&        srcTime,
                       const GObservation& obs,
                       const bool&         grad) const;
    GVector  eval_probs(const GModelSky&    model,
                        const GObservation& obs,
                        GMatrixSparse*      gradients = NULL) const;
    int      size_edisp_vector(const GModelSky&    model,
                               const GObservation& obs,
                               const bool&         grad) const;
    GEbounds ebounds_model(const GModelSky& model) const;

    // Virtual protected methods
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
    virtual double irf_composite(const GEvent&       event,
                                 const GSource&      source,
                                 const GObservation& obs) const;
    virtual GVector irf_ptsrc(const GModelSky&    model,
                              const GObservation& obs,
                              GMatrix*            gradients = NULL) const;
    virtual GVector irf_radial(const GModelSky&    model,
                               const GObservation& obs,
                               GMatrix*            gradients = NULL) const;
    virtual GVector irf_elliptical(const GModelSky&    model,
                                   const GObservation& obs,
                                   GMatrix*            gradients = NULL) const;
    virtual GVector irf_diffuse(const GModelSky&    model,
                                const GObservation& obs,
                                GMatrix*            gradients = NULL) const;
    virtual GVector irf_composite(const GModelSky&    model,
                                  const GObservation& obs,
                                  GMatrix*            gradients = NULL) const;

    // Protected classes
    class edisp_kerns : public GFunctions {
    public:
        edisp_kerns(const GResponse*    parent,
                    const GObservation* obs,
                    const GModelSky*    model,
                    const GEvent*       event,
                    const GTime&        srcTime,
                    const bool&         grad);
        int     size(void) const { return m_size; }
        GVector eval(const double& etrue);
    protected:
        const GResponse*        m_parent;  //!< Response
        const GObservation*     m_obs;     //!< Observation
        const GModelSky*        m_model;   //!< Sky model
        const GEvent*           m_event;   //!< Event
        int                     m_size;    //!< Array of values and gradients
        std::vector<GModelPar*> m_pars;    //!< Parameter pointers
        GTime                   m_srcTime; //!< True arrival time
        bool                    m_grad;    //!< Gradient flag
    };
    class irf_radial_kern_theta : public GFunction {
    public:
        irf_radial_kern_theta(const GResponse*           rsp,
                              const GEvent*              event,
                              const GObservation*        obs,
                              const GModelSpatialRadial* model,
                              const GMatrix*             rot,
                              const GEnergy*             srcEng,
                              const GTime*               srcTime,
                              int                        iter_phi) :
                              m_rsp(rsp),
                              m_event(event),
                              m_obs(obs),
                              m_model(model),
                              m_rot(rot),
                              m_srcEng(srcEng),
                              m_srcTime(srcTime),
                              m_iter_phi(iter_phi) { }
        double eval(const double& phi);
    protected:
        const GResponse*           m_rsp;       //!< Response
        const GEvent*              m_event;     //!< Event
        const GObservation*        m_obs;       //!< Observation
        const GModelSpatialRadial* m_model;     //!< Radial model
        const GMatrix*             m_rot;       //!< Rotation matrix
        const GEnergy*             m_srcEng;    //!< True photon energy
        const GTime*               m_srcTime;   //!< Arrival time
        int                        m_iter_phi;  //!< Iterations in phi
    };
    class irf_radial_kern_phi : public GFunction {
    public:
        irf_radial_kern_phi(const GResponse*    rsp,
                            const GEvent*       event,
                            const GObservation* obs,
                            const GMatrix*      rot,
                            const double&       theta,
                            const GEnergy*      srcEng,
                            const GTime*        srcTime) :
                            m_rsp(rsp),
                            m_event(event),
                            m_obs(obs),
                            m_rot(rot),
                            m_sin_theta(std::sin(theta)),
                            m_cos_theta(std::cos(theta)),
                            m_srcEng(srcEng),
                            m_srcTime(srcTime) { }
        double eval(const double& phi);
    protected:
        const GResponse*    m_rsp;       //!< Response
        const GEvent*       m_event;     //!< Event
        const GObservation* m_obs;       //!< Observation
        const GMatrix*      m_rot;       //!< Rotation matrix
        double              m_sin_theta; //!< sin(theta)
        double              m_cos_theta; //!< cos(theta)
        const GEnergy*      m_srcEng;    //!< True photon energy
        const GTime*        m_srcTime;   //!< Arrival time
    };
    class irf_elliptical_kern_theta : public GFunction {
    public:
        irf_elliptical_kern_theta(const GResponse*               rsp,
                                  const GEvent*                  event,
                                  const GObservation*            obs,
                                  const GModelSpatialElliptical* model,
                                  const GMatrix*                 rot,
                                  const GEnergy*                 srcEng,
                                  const GTime*                   srcTime,
                                  int                            iter_phi) :
                                  m_rsp(rsp),
                                  m_event(event),
                                  m_obs(obs),
                                  m_model(model),
                                  m_rot(rot),
                                  m_srcEng(srcEng),
                                  m_srcTime(srcTime),
                                  m_iter_phi(iter_phi) { }
        double eval(const double& phi);
    protected:
        const GResponse*               m_rsp;       //!< Response
        const GEvent*                  m_event;     //!< Event
        const GObservation*            m_obs;       //!< Observation
        const GModelSpatialElliptical* m_model;     //!< Elliptical model
        const GMatrix*                 m_rot;       //!< Rotation matrix
        const GEnergy*                 m_srcEng;    //!< True photon energy
        const GTime*                   m_srcTime;   //!< Arrival time
        int                            m_iter_phi;  //!< Iterations in phi
    };
    class irf_elliptical_kern_phi : public GFunction {
    public:
        irf_elliptical_kern_phi(const GResponse*               rsp,
                                const GEvent*                  event,
                                const GObservation*            obs,
                                const GModelSpatialElliptical* model,
                                const GMatrix*                 rot,
                                const double&                  theta,
                                const GEnergy*                 srcEng,
                                const GTime*                   srcTime) :
                                m_rsp(rsp),
                                m_event(event),
                                m_obs(obs),
                                m_model(model),
                                m_rot(rot),
                                m_theta(theta),
                                m_sin_theta(std::sin(theta)),
                                m_cos_theta(std::cos(theta)),
                                m_srcEng(srcEng),
                                m_srcTime(srcTime) { }
        double eval(const double& phi);
    protected:
        const GResponse*               m_rsp;       //!< Response
        const GEvent*                  m_event;     //!< Event
        const GObservation*            m_obs;       //!< Observation
        const GModelSpatialElliptical* m_model;     //!< Elliptical model
        const GMatrix*                 m_rot;       //!< Rotation matrix
        double                         m_theta;     //!< Theta
        double                         m_sin_theta; //!< sin(theta)
        double                         m_cos_theta; //!< cos(theta)
        const GEnergy*                 m_srcEng;    //!< True photon energy
        const GTime*                   m_srcTime;   //!< Arrival time
    };

    // Protected members
    bool   m_use_irf_cache;             //!< Control usage of irf cache
    bool   m_use_nroi_cache;            //!< Control usage of nroi cache
    int    m_irf_radial_iter_theta;     //!< Radial model integration theta iterations
    int    m_irf_radial_iter_phi;       //!< Radial model integration phi iterations
    int    m_irf_elliptical_iter_theta; //!< Elliptical model integration theta iterations
    int    m_irf_elliptical_iter_phi;   //!< Elliptical model integration phi iterations
    double m_irf_diffuse_resolution;    //!< Angular resolution for diffuse model

    // Cache for irf(GEvent&, GSource&, GObservation&) and
    // nroi(GModelSky&, GEnergy&, GTime&, GEnergy&, GTime&, GObservation&)
    // computations
    mutable GResponseCache       m_irf_cache;
    mutable GResponseCache       m_nroi_cache;
    mutable GResponseVectorCache m_irf_vector_cache;
};

#endif /* GRESPONSE_HPP */
