/***************************************************************************
 *                  GEvent.hpp  -  Event abstract base class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEvent.hpp
 * @brief GEvent abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENT_HPP
#define GEVENT_HPP

/* __ Includes ___________________________________________________________ */
#include "GModels.hpp"
#include "GVector.hpp"
#include "GInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPointing.hpp"
#include "GResponse.hpp"


/***********************************************************************//**
 * @class GEvent
 *
 * @brief Abstract interface for the event classes.
 *
 * This class provides an abstract interface to a event. A event can be
 * either a physical event occuring in the instrument (called an event
 * atom) or a collection of events with similar properties (called an
 * event bin). While event atoms are used for unbinned analysis, event bins
 * are used for binned analysis. The methods isatom() and isbin() inform
 * whether an event is either an atom or a bin.
 * The counts() method returns the number of event atoms in a event bin.
 * For an event atom, this method returns by definition 1.
 * The model() method returns either the probability for an event atom
 * to occur (unbinned analysis) or the expected number of counts for an
 * event bin (binned analysis).
 * Attributes of an event atom or bin can be accessed through the dir(),
 * energy(), time(), pnt(), rsp() methods that all return const pointers
 * to the relevant information.
 *
 * This method does not hold any data members. Data members are stored in
 * the derived classes.
 ***************************************************************************/
class GEvent {

    // Friend classes
    friend class GEvents;

public:
    // Constructors and destructors
    GEvent(void);
    GEvent(const GEvent& event);
    virtual ~GEvent(void);

    // Operators
    virtual GEvent& operator= (const GEvent& event);

    // Virtual methods
    virtual bool             isatom(void) const = 0;
    virtual bool             isbin(void) const = 0;
    virtual double           counts(void) const = 0;
    virtual double           model(GModels& models, GVector* gradient = NULL) const = 0;
    virtual const GInstDir*  dir(void) const = 0;
    virtual const GEnergy*   energy(void) const = 0;
    virtual const GTime*     time(void) const = 0;
    virtual const GPointing* pnt(void) const = 0;
    virtual const GResponse* rsp(void) const = 0;
    
protected:
    // Protected methods
    void            init_members(void);
    void            copy_members(const GEvent& event);
    void            free_members(void);
    virtual GEvent* clone(void) const = 0;

private:
};

#endif /* GEVENT_HPP */
