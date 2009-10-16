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


/***********************************************************************//**
 * @class GEvent
 *
 * @brief Abstract interface for the event classes.
 ***************************************************************************/
class GEvent {

    // Friend classes
    friend class GEvents;

public:
    // Constructors and destructors
    GEvent();
    GEvent(const GEvent& event);
    virtual ~GEvent();

    // Operators
    virtual GEvent& operator= (const GEvent& event);

    // Virtual methods
    virtual double counts(void) const = 0;
    virtual double model(GModels& models, GVector* gradient) const = 0;
    virtual bool   isatom(void) const = 0;
    virtual bool   isbin(void) const = 0;
    
    // Implemented methods
    
protected:
    // Protected methods
    void            init_members(void);
    void            copy_members(const GEvent& event);
    void            free_members(void);
    virtual GEvent* clone(void) const = 0;

    // Protected data area

private:
};

#endif /* GEVENT_HPP */
