/***************************************************************************
 *           GEventList.hpp  -  Abstract event list container class        *
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
 * @file GEventList.hpp
 * @brief GEventList container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTLIST_HPP
#define GEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvents.hpp"
#include "GEventAtom.hpp"
#include "GFits.hpp"


/***********************************************************************//**
 * @class GEventList
 *
 * @brief GEventList container class interface defintion.
 ***************************************************************************/
class GEventList : public GEvents {

	// Friend classes
    friend class GData;
    friend class GObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEventList& list);

public:
    // Constructors and destructors
    GEventList();
    GEventList(const GEventList& list);
    virtual ~GEventList();

    // Operators
    virtual GEventList& operator= (const GEventList& list);

    // Virtual methods
	virtual void        load(const std::string& filename) = 0;
    virtual GEventAtom* pointer(int index) = 0;

    // Implemented methods
    int  number(void) const;
    int  elements(void) const;
    bool islist(void) const { return true; }
    bool iscube(void) const { return false; }

protected:
    // Protected methods
    void                init_members(void);
    void                copy_members(const GEventList& list);
    void                free_members(void);
    virtual GEventList* clone(void) const = 0;

    // Protected data area
    int         m_num;           //!< Number of events
    GEventAtom* m_events;        //!< Pointer to events

private:
};

#endif /* GEVENTLIST_HPP */
