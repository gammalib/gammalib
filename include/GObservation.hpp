/***************************************************************************
 *           GObservation.hpp  -  Observation abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservation.hpp
 * @brief GObservation abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOBSERVATION_HPP
#define GOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvents.hpp"
#include "GResponse.hpp"
#include "GRoi.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GModels.hpp"
#include "GTime.hpp"
#include "GEnergy.hpp"
#include "GIntegrand.hpp"

/* __ Typedefs ___________________________________________________________ */


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract interface for the observation classes.
 *
 * The class implements support methods to compute Npred values and
 * parameter gradients. These support methods can be overloaded by derived
 * classes that implement instrument specific observations in order to
 * optimize the execution speed for model fitting.
 ***************************************************************************/
class GObservation {

    // Friend classes
    friend class GObservations;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GObservation& obs);

public:
    // Constructors and destructors
    GObservation(void);
    GObservation(const GObservation& obs);
    virtual ~GObservation(void);

    // Operators
    virtual GObservation& operator= (const GObservation& obs);

    // Pure virtual methods
    virtual void   response(const std::string& irfname, std::string caldb = "") = 0;

    // Virtual methods
    virtual double model(const GModels& models, const GInstDir& obsDir,
                         const GEnergy& obsEng, const GTime& obsTime,
                         const GPointing& pnt, GVector* gradient = NULL) const;
    virtual double npred(const GModels& models, GVector* gradient = NULL) const;

    // Implemented methods
    void        obsname(const std::string& obsname) { m_obsname=obsname; return; }
    void        instrument(const std::string& instrument) { m_instrument=instrument; return; }
    void        roi(const GRoi& roi) { m_roi=roi.clone(); return; }
    void        ebounds(const GEbounds& ebounds) { m_ebounds=ebounds; return; }
    void        gti(const GGti& gti) { m_gti=gti; return; }
    std::string obsname(void) const { return m_obsname; }
    std::string instrument(void) const { return m_instrument; }
    GTime       tstart(void) const { return m_gti.tstart(); }
    GTime       tstop(void) const { return  m_gti.tstop(); }
    GEnergy     emin(void) const { return m_ebounds.emin(); }
    GEnergy     emax(void) const { return m_ebounds.emax(); }
    GRoi*       roi(void) { return m_roi; }
    GEbounds*   ebounds(void) { return &m_ebounds; }
    GGti*       gti(void) { return &m_gti; }
    GEvents*    events(void) const { return m_events; }
    GResponse*  response(void) const { return m_response; }

protected:
    // Protected methods
    void                  init_members(void);
    void                  copy_members(const GObservation& obs);
    void                  free_members(void);
    virtual GObservation* clone(void) const = 0;

    // Npred integration methods
    virtual double npred_kern(const GModel& model, const GSkyDir& srcDir,
                              const GEnergy& srcEng, const GTime& srcTime,
                              const GPointing& pnt) const;
    virtual double npred_spat(const GModel& model,
                              const GEnergy& srcEng, const GTime& srcTime,
                              const GPointing& pnt) const;
    virtual double npred_spec(const GModel& model, const GTime& srcTime) const;
    virtual double npred_temp(const GModel& model) const;

    // Npred gradient integration methods
    virtual double npred_grad_kern(const GModel& model, int ipar,
                                   const GSkyDir& srcDir, const GEnergy& srcEng,
                                   const GTime& srcTime, const GPointing& pnt) const;
    virtual double npred_grad_spat(const GModel& model, int ipar,
                                   const GEnergy& srcEng, const GTime& srcTime,
                                   const GPointing& pnt) const;
    virtual double npred_grad_spec(const GModel& model, int ipar,
                                   const GTime& srcTime) const;
    virtual double npred_grad_temp(const GModel& model, int ipar) const;

    // Npred kernel classes
    class npred_kern_spat : public GIntegrand {
    public:
        npred_kern_spat(const GObservation* parent, const GModel& model,
                        const GTime& srcTime, const GPointing* pnt) :
                        m_parent(parent), m_model(&model), m_time(&srcTime),
                        m_pnt(pnt) { return; }
        double eval(double x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to model
        const GTime*        m_time;   //!< Pointer to time
        const GPointing*    m_pnt;    //!< Pointer to telescope pointing
    };
    class npred_kern_spec : public GIntegrand {
    public:
        npred_kern_spec(const GObservation* parent, const GModel& model) :
                        m_parent(parent), m_model(&model) { return; }
        double eval(double x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to model
    };

    // Npred gradient kernel classes
    class npred_grad_kern_spat : public GIntegrand {
    public:
        npred_grad_kern_spat(const GObservation* parent, const GModel& model,
                             int ipar, const GTime& srcTime,
                             const GPointing* pnt) :
                             m_parent(parent), m_model(&model), m_ipar(ipar),
                             m_time(&srcTime), m_pnt(pnt) { return; }
        double eval(double x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to model
        int                 m_ipar;   //!< Parameter index
        const GTime*        m_time;   //!< Pointer to time
        const GPointing*    m_pnt;    //!< Pointer to pointing
    };
    class npred_grad_kern_spec : public GIntegrand {
    public:
        npred_grad_kern_spec(const GObservation* parent, const GModel& model,
                             int ipar) : m_parent(parent), m_model(&model),
                             m_ipar(ipar) { return; }
        double eval(double x);
    protected:
        const GObservation* m_parent; //!< Pointer to parent
        const GModel*       m_model;  //!< Pointer to model
        int                 m_ipar;   //!< Parameter index
    };

    // Protected data area
    std::string m_obsname;      //!< Name of observation
    std::string m_instrument;   //!< Instrument name
    GEbounds    m_ebounds;      //!< Energy boundaries used for analysis
    GGti        m_gti;          //!< Good time intervals used for analysis
    GRoi*       m_roi;          //!< Pointer to region of interest used for analysis
    GEvents*    m_events;       //!< Pointer to events
    GResponse*  m_response;     //!< Pointer to instrument response functions

private:
};

#endif /* GOBSERVATION_HPP */
