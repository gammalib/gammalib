/***************************************************************************
 *                GCTAEventList.hpp  -  CTA Event list class               *
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
 * @file GCTAEventList.hpp
 * @brief GCTAEventList class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAEVENTLIST_HPP
#define GCTAEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventList.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAObservation.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GCTAEventList
 *
 * @brief GCTAEventList class interface defintion.
 ***************************************************************************/
class GCTAEventList : public GEventList {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAEventList& list);

public:
    // Constructors and destructors
    GCTAEventList(void);
    GCTAEventList(const GCTAEventList& list);
    virtual ~GCTAEventList(void);

    // Operators
    GCTAEventList& operator= (const GCTAEventList& list);

    // Implemented pure virtual methods
    void           clear(void);
    GCTAEventList* clone(void) const;
    void           load(const std::string& filename);
    GCTAEventAtom* pointer(int index);
    int            number(void) const { return m_num; }
    int            size(void) const { return m_num; }

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GCTAEventList& list);
    void           free_members(void);
    void           load_events(GFitsTable* hdu);

    // Protected data area
    int            m_num;            //!< Number of events
    GCTAEventAtom* m_events;         //!< Pointer to events
};

#endif /* GCTAEVENTLIST_HPP */
