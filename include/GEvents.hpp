/***************************************************************************
 *                 GEvents.hpp  -  Events container class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEvents.hpp
 * @brief GEvents container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTS_HPP
#define GEVENTS_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvent.hpp"
#include "GFits.hpp"


/***********************************************************************//**
 * @class GEvents
 *
 * @brief GEvents container class interface defintion.
 *
 * This class is an abstract container base class for events. Events are
 * generally associated to an observation, and the class keeps the pointer
 * to an existing observations as a member.
 ***************************************************************************/
class GEvents {

	// Friend classes
    friend class GObservations;
    friend class GObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEvents& events);

public:
    // Constructors and destructors
    GEvents();
    GEvents(const GEvents& events);
    virtual ~GEvents();

    // Operators
    virtual GEvents& operator= (const GEvents& events);

    // Virtual methods
	virtual void    load(const std::string& filename) = 0;
    virtual GEvent* pointer(int index) = 0;
    virtual int     number(void) const = 0;
    virtual int     size(void) const = 0;
    virtual bool    islist(void) const = 0;
    virtual bool    iscube(void) const = 0;
    
    // Other methods
    void          obs(GObservation* ptr) { m_obs=ptr; return; }
    GObservation* obs(void) const { return m_obs; }
    
    // Event iterator
    class iterator {
    friend class GEvents;
    public:
        iterator();
        iterator(GEvents *events);
        ~iterator() { return; }
        iterator& operator++(void) { m_index++; return *this; }
        iterator  operator++(int junk);
        bool      operator==(const iterator& it) const { return (m_index == it.m_index); }
        bool      operator!=(const iterator& it) const { return (m_index != it.m_index); }
        GEvent&   operator*(void) { return *(m_base->pointer(m_index)); }
        GEvent*   operator->(void) { return m_base->pointer(m_index); }
    protected:
        int      m_index;        //!< Actuel event index
        int      m_num;          //!< Number of events in GEvents object
        GEvents *m_base;         //!< Pointer to GEvents object
    };
    iterator begin(void);
    iterator end(void);

protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GEvents& events);
    void             free_members(void);
    virtual GEvents* clone(void) const = 0;
    
    // Protected data
    GObservation*  m_obs;        //!< Pointer to associated observation
    
};

#endif /* GEVENTS_HPP */
