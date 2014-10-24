/***************************************************************************
 *            GObservation.hpp - Abstract observation base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
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
 * @file GObservation.hpp
 * @brief Abstract observation base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GOBSERVATION_HPP
#define GOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GEvents.hpp"
#include "GResponse.hpp"
#include "GModels.hpp"
#include "GTime.hpp"
#include "GEnergy.hpp"
#include "GFunction.hpp"
#include "GVector.hpp"
#include "GMatrixSparse.hpp"


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract observation base class
 *
 * This class provides an abstract interface for an observation. The
 * observation collects information about the instrument, holds the measured
 * events, and provides information about the analysis definition.
 *
 * The response() method returns a pointer to the response function. The
 * derived classes have to make sure that this method never returns NULL.
 *
 * The method model() returns the probability for an event to be measured
 * with a given instrument direction, a given energy and at a given time,
 * given a source model and an instrument pointing direction.
 * The method npred() returns the total number of expected events within the
 * analysis region for a given source model and a given instrument pointing
 * direction.
 * The methods are defined as virtual and can be overloaded by derived classes
 * that implement instrument specific observations in order to optimize the
 * execution speed for data analysis.
 ***************************************************************************/
class GObservation : public GBase {

public:
    // Constructors and destructors
    GObservation(void);
    GObservation(const GObservation& obs);
    virtual ~GObservation(void);

    // Operators
    virtual GObservation& operator=(const GObservation& obs);

    // Pure virtual methods
    virtual void             clear(void) = 0;
    virtual GObservation*    clone(void) const = 0;
    virtual std::string      classname(void) const = 0;
    virtual void             response(const GResponse& rsp) = 0;
    virtual const GResponse* response(void) const = 0;
    virtual std::string      instrument(void) const = 0;
    virtual double           ontime(void) const = 0;
    virtual double           livetime(void) const = 0;
    virtual double           deadc(const GTime& time) const = 0;
    virtual void             read(const GXmlElement& xml) = 0;
    virtual void             write(GXmlElement& xml) const = 0;
    virtual std::string      print(const GChatter& chatter = NORMAL) const = 0;

    // Virtual methods
    virtual const GEvents*   events(void) const;
    virtual void             events(const GEvents& events);
    virtual double           likelihood(const GModels& models,
                                        GVector*       gradient,
                                        GMatrixSparse* curvature,
                                        double*        npred) const;
    virtual double           model(const GModels& models,
                                   const GEvent&  event,
                                   GVector*       gradient = NULL) const;
    virtual double           npred(const GModels& models,
                                   GVector*       gradient = NULL) const;
    virtual double           npred(const GModel& model) const;
    virtual double           model_grad(const GModel&    model,
                                        const GModelPar& par,
                                        const GEvent&    event) const;
    virtual double           npred_grad(const GModel&    model,
                                        const GModelPar& par) const;

    // Implemented methods
    void               name(const std::string& name);
    void               id(const std::string& id);
    void               statistics(const std::string& statistics);
    const std::string& name(void) const;
    const std::string& id(void) const;
    const std::string& statistics(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GObservation& obs);
    void free_members(void);

    // Likelihood methods
    virtual double likelihood_poisson_unbinned(const GModels& models,
                                               GVector*       gradient,
                                               GMatrixSparse* curvature,
                                               double*        npred) const;
    virtual double likelihood_poisson_binned(const GModels& models,
                                             GVector*       gradient,
                                             GMatrixSparse* curvature,
                                             double*        npred) const;
    virtual double likelihood_gaussian_binned(const GModels& models,
                                              GVector*       gradient,
                                              GMatrixSparse* curvature,
                                              double*        npred) const;

    // Model gradient kernel classes
    class model_func : public GFunction {
    public:
        model_func(const GObservation* parent,
                   const GModel&       model,
                   GModelPar*          par,
                   const GEvent&       event) :
                   m_parent(parent),
                   m_model(model),
                   m_par(par),
                   m_event(event) { }
        double eval(const double& x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel&       m_model;  //!< Reference to model
        GModelPar*          m_par;    //!< Pointer to parameter
        const GEvent&       m_event;  //!< Reference to event
    };

    // Npred methods
    virtual double npred_spec(const GModel& model, const GTime& obsTime) const;

    // Npred kernel classes
    class npred_kern : public GFunction {
    public:
        npred_kern(const GObservation* parent,
                   const GModel*       model) :
                   m_parent(parent),
                   m_model(model) { }
        double eval(const double& x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to source model
    };

    class npred_spec_kern : public GFunction {
    public:
        npred_spec_kern(const GObservation* parent,
                        const GModel*       model,
                        const GTime*        obsTime) :
                        m_parent(parent),
                        m_model(model),
                        m_time(obsTime) { }
        double eval(const double& x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to source model
        const GTime*        m_time;   //!< Pointer to time
    };

    // Npred gradient kernel classes
    class npred_func : public GFunction {
    public:
        npred_func(const GObservation* parent,
                   const GModel&       model,
                   GModelPar*          par) :
                   m_parent(parent),
                   m_model(model),
                   m_par(par) { }
        double eval(const double& x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel&       m_model;  //!< Reference to model
        GModelPar*          m_par;    //!< Pointer to parameter
    };

    // Protected data area
    std::string m_name;        //!< Observation name
    std::string m_id;          //!< Observation identifier
    std::string m_statistics;  //!< Optimizer statistics (default=Poisson)
    GEvents*    m_events;      //!< Pointer to event container
};


/***********************************************************************//**
 * @brief Set observation name
 *
 * @param[in] name Observation name.
 *
 * Set name of the observation.
 ***************************************************************************/
inline
void GObservation::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Set observation identifier
 *
 * @param[in] id Observation identifier.
 *
 * Set identifier of the observation.
 ***************************************************************************/
inline
void GObservation::id(const std::string& id)
{
    m_id = id;
    return;
}


/***********************************************************************//**
 * @brief Set optimizer statistics
 *
 * @param[in] statistics Optimizer statistics.
 *
 * Set optimizer statistics for the observation.
 ***************************************************************************/
inline
void GObservation::statistics(const std::string& statistics)
{
    m_statistics = statistics;
    return;
}


/***********************************************************************//**
 * @brief Return observation name
 *
 * @return Observation name.
 ***************************************************************************/
inline
const std::string& GObservation::name(void) const
{
    return (m_name);
}


/***********************************************************************//**
 * @brief Return observation identifier
 *
 * @return Observation identifier.
 ***************************************************************************/
inline
const std::string& GObservation::id(void) const
{
    return (m_id);
}


/***********************************************************************//**
 * @brief Return optimizer statistics
 *
 * @return Optimizer statistics.
 ***************************************************************************/
inline
const std::string& GObservation::statistics(void) const
{
    return (m_statistics);
}

#endif /* GOBSERVATION_HPP */
