/***************************************************************************
 *           GObservation.hpp  -  Observation abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2009 by Jurgen Knodlseder                           *
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
#include "GGti.hpp"
#include "GModels.hpp"

/* __ Typedefs ___________________________________________________________ */


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract interface for the observation classes.
 ***************************************************************************/
class GObservation {

    // Friend classes
    friend class GData;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GObservation& obs);

public:
    // Constructors and destructors
    GObservation();
    GObservation(const GObservation& obs);
    virtual ~GObservation();

    // Operators
    virtual GObservation& operator= (const GObservation& obs);

    // Methods
    void         obsname(const std::string& obsname) { m_obsname=obsname; return; }
    void         instrument(const std::string& instrument) { m_instrument=instrument; return; }
    void         tstart(const double& tstart) { m_tstart=tstart; return; }
    void         tstop(const double& tstop) { m_tstop=tstop; return; }
    void         emin(const double& emin) { m_emin=emin; return; }
    void         emax(const double& emax) { m_emax=emax; return; }
    void         gti(const GGti& gti) { m_gti=gti; return; }
    double       tstart(void) const { return m_tstart; }
    double       tstop(void) const { return m_tstop; }
    double       emin(void) const { return m_emin; }
    double       emax(void) const { return m_emax; }
    std::string  obsname(void) const { return m_obsname; }
    std::string  instrument(void) const { return m_instrument; }
    GEvents*     events(void) const { return m_events; }
    GResponse*   response(void) const { return m_response; }
    GGti*        gti(void) { return &m_gti; }

protected:
    // Protected methods
    void                  init_members(void);
    void                  copy_members(const GObservation& obs);
    void                  free_members(void);
    virtual GObservation* clone(void) const = 0;

    // Protected data area
    std::string m_obsname;      //!< Name of observation
    std::string m_instrument;   //!< Instrument name
    double      m_tstart;       //!< Start time of observation
    double      m_tstop;        //!< Stop time of observations
    double      m_emin;         //!< Minimum energy covered by observation
    double      m_emax;         //!< Maximum energy covered by observation
    GEvents*    m_events;       //!< Pointer to events
    GResponse*  m_response;     //!< Pointer to instrument response functions
    GGti        m_gti;          //!< Good time intervals

private:
};

#endif /* GOBSERVATION_HPP */
