/***************************************************************************
 *           GObservation.hpp  -  Abstract observation base class          *
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
#include "GPointing.hpp"
#include "GModels.hpp"
#include "GTime.hpp"
#include "GEnergy.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract interface for the observation classes
 *
 * This class provides an abstract interface for an observation. The
 * observation collects information about the instrument, holds the measured
 * events, and provides information about the analysis definiton.
 * The method model() returns the probability for an event to be measured
 * with a given instrument direction, a given energy and at a given time,
 * given a source model and an instrument pointing direction.
 * The method npred() returns the total number of expected events within the
 * analysis region for a given source model and a given instrument pointing
 * direction.
 * The methods a defined as virtual and can be overloaded by derived classes
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
    virtual GObservation& operator= (const GObservation& obs);

    // Pure virtual methods
    virtual void          clear(void) = 0;
    virtual GObservation* clone(void) const = 0;
    virtual void          response(const GResponse& rsp) = 0;
    virtual GResponse*    response(void) const = 0;
    virtual GPointing*    pointing(void) const = 0;
    virtual std::string   instrument(void) const = 0;
    virtual double        ontime(void) const = 0;
    virtual double        livetime(void) const = 0;
    virtual double        deadc(const GTime& time) const = 0;
    virtual void          read(const GXmlElement& xml) = 0;
    virtual void          write(GXmlElement& xml) const = 0;
    virtual std::string   print(void) const = 0;

    // Virtual methods
    virtual double        model(const GModels& models, const GEvent& event,
                                GVector* gradient = NULL) const;
    virtual double        npred(const GModels& models, GVector* gradient = NULL) const;

    // Implemented methods
    void                  name(const std::string& name);
    void                  id(const std::string& id);
    void                  events(const GEvents* events);
    void                  statistics(const std::string& statistics);
    const std::string&    name(void) const { return m_name; }
    const std::string&    id(void) const { return m_id; }
    const GEvents*        events(void) const;
    const std::string&    statistics(void) const { return m_statistics; }

    // Other methods
    virtual double model_grad(const GModel& model, const GEvent& event, int ipar) const;
    virtual double npred_grad(const GModel& model, int ipar) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GObservation& obs);
    void free_members(void);


    // Model gradient kernel classes
    class model_func : public GFunction {
    public:
        model_func(const GObservation* parent,
                   const GModel&       model,
                   const GEvent&       event,
                   int                 ipar) :
                   m_parent(parent),
                   m_model(&model),
                   m_event(&event),
                   m_ipar(ipar) { }
        double eval(double x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to model
        const GEvent*       m_event;  //!< Pointer to event
        int                 m_ipar;   //!< Parameter index
    };

    // Npred methods
    virtual double npred_temp(const GModel& model) const;
    virtual double npred_spec(const GModel& model, const GTime& obsTime) const;

    // Npred kernel classes
    class npred_temp_kern : public GFunction {
    public:
        npred_temp_kern(const GObservation* parent,
                        const GModel*       model) :
                        m_parent(parent),
                        m_model(model) { }
        double eval(double x);
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
        double eval(double x);
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
                   int                 ipar) :
                   m_parent(parent),
                   m_model(&model),
                   m_ipar(ipar) { return; }
        double eval(double x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to model
        int                 m_ipar;   //!< Parameter index
    };

    // Protected data area
    std::string m_name;         //!< Name of observation
    std::string m_id;           //!< Observation identifier
    std::string m_statistics;   //!< Optimizer statistics (default=poisson)
    GEvents*    m_events;       //!< Pointer to event container
};

#endif /* GOBSERVATION_HPP */
