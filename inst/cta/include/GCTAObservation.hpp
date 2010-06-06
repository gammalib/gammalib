/***************************************************************************
 *               GCTAObservation.hpp  -  CTA Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAObservation.hpp
 * @brief GCTAObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAOBSERVATION_HPP
#define GCTAOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GIntegrand.hpp"


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief Interface class for CTA observations.
 ***************************************************************************/
class GCTAObservation : public GObservation {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAObservation& obs);

public:
    // Constructors and destructors
    GCTAObservation(void);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Operators
    GCTAObservation& operator= (const GCTAObservation& obs);

    // Methods
    void   response(const std::string& irfname, std::string caldb = "");
    void   load_unbinned(const std::string& evname);
    double npred(const GModels& pars) const;

protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GCTAObservation& obs);
    void             free_members(void);
    GCTAObservation* clone(void) const;
    
    //
    double npred_integrand(const GModel& model, const GSkyDir& srcDir,
                           const GEnergy& srcEng, const GTime& srcTime,
                           const GPointing& pnt) const;
    double npred_integrate_spatial(const GModel& model,
                                   const GEnergy& srcEng,
                                   const GTime& srcTime,
                                   const GPointing& pnt) const;
    double npred_integrate_spectral(const GModel& model,
                                    const GTime& srcTime,
                                    const GPointing& pnt) const;
    double npred_integrate_temporal(const GModel& model,
                                    const GPointing& pnt) const;

    // Spectral integrand
    class int_spec : public GIntegrand {
    public:
        int_spec(const GCTAObservation* parent, const GModel& model, 
                 const GTime& srcTime, const GPointing& pnt) : 
                 m_parent(parent), m_model(&model), m_time(&srcTime), 
                 m_pnt(&pnt) { return; }
        double eval(double x) { 
                 GEnergy eng;
                 eng.TeV(x);
                 return (m_parent->npred_integrate_spatial(*m_model,eng,*m_time,*m_pnt));
                 }
    protected:
        const GCTAObservation* m_parent; //!< Pointer to parent
        const GModel*          m_model;  //!< Pointer to model
        const GTime*           m_time;   //!< Pointer to time
        const GPointing*       m_pnt;    //!< Pointer to pointing
    };

    // Temporal integrand
    class int_temp : public GIntegrand {
    public:
        int_temp(const GCTAObservation* parent, const GModel& model, const GPointing& pnt) : 
                 m_parent(parent), m_model(&model), m_pnt(&pnt) { return; }
        double eval(double x) { 
                 GTime time;
                 time.met(x);
                 return (m_parent->npred_integrate_spectral(*m_model,time,*m_pnt));
                 }
    protected:
        const GCTAObservation* m_parent; //!< Pointer to parent
        const GModel*          m_model;  //!< Pointer to model
        const GPointing*       m_pnt;    //!< Pointer to pointing
    };
    
    // Protected data area
};

#endif /* GCTAOBSERVATION_HPP */
