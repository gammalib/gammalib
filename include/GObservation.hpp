/***************************************************************************
 *           GObservation.hpp  -  Observation abstract base class          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
#include "GResponse.hpp"
#include "GEvents.hpp"
#include "GGti.hpp"


/* __ Typedefs ___________________________________________________________ */
typedef void (GObservation::*fct_free)(void);


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
	void         link(void);
	void         release(void);
	double       tstart(void);
	double       tstop(void);
	std::string  obsname(void);
	std::string  instrument(void);
  
protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GObservation& obs);
    void    free_members(void);
    virtual GObservation* clone(void) const = 0;

    // Protected data area
    std::string m_obsname;      //!< Name of observation
	std::string m_instrument;   //!< Instrument name
	double      m_tstart;       //!< Start time of observation
	double      m_tstop;        //!< Stop time of observations
	GEvents*    m_events;       //!< Pointer to events
	GGti*       m_gti;          //!< Pointer to good time intervals
	GResponse*  m_response;     //!< Pointer to instrument response functions
	int         m_num_links;    //!< Number of links into GData objects
	int         m_free_pending; //!< Signals pending destruction (1=pending, 0=no pending)
	fct_free    m_free_members; //!< Pointer to derive class free member function

private:
};

#endif /* GOBSERVATION_HPP */
