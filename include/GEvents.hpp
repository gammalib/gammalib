/***************************************************************************
 *                 GEvents.hpp  -  Events container class                  *
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
 ***************************************************************************/
class GEvents {

	// Friend classes
    friend class GData;
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
//    virtual void    load(GFitsHDU* hdu) = 0;
    virtual GEvent* pointer(int index) const = 0;
    virtual int     number(void) const = 0;
    virtual int     elements(void) const = 0;
    
    // Implemented methods

    // Event iterator
    class iterator {
    friend class GEvents;
    public:
        iterator();
        iterator(GEvents *events);
        ~iterator();
        iterator& operator++(void);                // Prefix
        iterator  operator++(int junk);            // Postfix
        bool      operator==(const iterator& it) const;
        bool      operator!=(const iterator& it) const;
        GEvent&   operator*(void);
        GEvent*   operator->(void);
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
};

#endif /* GEVENTS_HPP */
