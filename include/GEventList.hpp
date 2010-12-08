/***************************************************************************
 *           GEventList.hpp  -  Abstract event list container class        *
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
 * @file GEventList.hpp
 * @brief GEventList container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTLIST_HPP
#define GEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEvents.hpp"
#include "GEventAtom.hpp"


/***********************************************************************//**
 * @class GEventList
 *
 * @brief GEventList container class interface defintion.
 ***************************************************************************/
class GEventList : public GEvents {

public:
    // Constructors and destructors
    GEventList(void);
    GEventList(const GEventList& list);
    virtual ~GEventList(void);

    // Operators
    virtual GEventList& operator= (const GEventList& list);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventList* clone(void) const = 0;
    virtual int         size(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual GEventAtom* pointer(int index) = 0;
    virtual int         number(void) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented pure virtual base class methods
    bool islist(void) const { return true; }
    bool iscube(void) const { return false; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEventList& list);
    void free_members(void);
};

#endif /* GEVENTLIST_HPP */
