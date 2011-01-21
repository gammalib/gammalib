/***************************************************************************
 *                GCTAEventList.hpp  -  CTA Event list class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
#include <string>
#include <vector>
#include "GEventList.hpp"
#include "GCTAEventAtom.hpp"
#include "GCTAObservation.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"


/***********************************************************************//**
 * @class GCTAEventList
 *
 * @brief GCTAEventList class interface defintion.
 ***************************************************************************/
class GCTAEventList : public GEventList {

public:
    // Constructors and destructors
    GCTAEventList(void);
    GCTAEventList(const GCTAEventList& list);
    virtual ~GCTAEventList(void);

    // Operators
    GCTAEventList& operator= (const GCTAEventList& list);

    // Implemented pure virtual base class methods
    void           clear(void);
    GCTAEventList* clone(void) const;
    int            size(void) const { return m_events.size(); }
    void           load(const std::string& filename);
    void           save(const std::string& filename, bool clobber = false) const;
    void           read(GFitsTable* hdu);
    void           write(GFits* file) const;
    GCTAEventAtom* pointer(int index);
    int            number(void) const { return m_events.size(); }
    std::string    print(void) const;

    // Implement other methods
    void append(const GCTAEventAtom& event);
    void reserve(const int& number);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAEventList& list);
    void free_members(void);
    void write_header(GFitsBinTable* hdu) const;

    // Protected members
    std::vector<GCTAEventAtom> m_events;  //!< Events
};

#endif /* GCTAEVENTLIST_HPP */
