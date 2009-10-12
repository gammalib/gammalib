/***************************************************************************
 *                  GEvent.hpp  -  Event abstract base class               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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


/***********************************************************************//**
 * @class GEvent
 *
 * @brief Abstract interface for the event classes.
 ***************************************************************************/
class GEvent {

    // Friend classes
    friend class GEvents;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEvent& event);

public:
    // Constructors and destructors
    GEvent();
    GEvent(const GEvent& event);
    virtual ~GEvent();

    // Operators
    virtual GEvent& operator= (const GEvent& event);

    // Virtual methods
    virtual std::string string(void) const = 0;
    
    // Implemented methods
    int test(void) const { return 47; }
    
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
