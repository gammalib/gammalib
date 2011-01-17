/***************************************************************************
 *     GEventList.i  -  Abstract event list container class python I/F     *
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
 * @file GEventList.i
 * @brief GEventList class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEventList.hpp"
%}


/***********************************************************************//**
 * @class GEventList
 *
 * @brief GEventList container class interface defintion.
 ***************************************************************************/
class GEventList : public GEvents {
public:
    // Constructors and destructors
    GEventList();
    GEventList(const GEventList& list);
    virtual ~GEventList();

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventList* clone(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual GEventAtom* pointer(int index) = 0;
    virtual int         number(void) const = 0;
    virtual int         size(void) const = 0;

    // Implemented pure virtul methods
    bool islist(void) const { return true; }
    bool iscube(void) const { return false; }
};


/***********************************************************************//**
 * @brief GEventList class extension
 *
 * The __getitem__ method makes the event list iteratable.
 ***************************************************************************/
%extend GEventList {
    GEventAtom* __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return self->pointer(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
};


/***********************************************************************//**
 * @brief GEventList type casts
 ***************************************************************************/
%inline %{
    GEventList* cast_GEventList(GEvents* events) {
        if (!events->islist())
            throw GException::fits_invalid_type("cast_GEventList(GEvents*)",
                                                "GEvents is not an event list.");
        return dynamic_cast<GEventList*>(events);
    }
%};
