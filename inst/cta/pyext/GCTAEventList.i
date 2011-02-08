/***************************************************************************
 *            GCTAEventList.i  -  CTA event atom container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAEventList.i
 * @brief CTA event atom container class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventList.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventList
 *
 * @brief CTA event atom container class Python interface
 ***************************************************************************/
class GCTAEventList : public GEventList {

public:
    // Constructors and destructors
    GCTAEventList(void);
    GCTAEventList(const GCTAEventList& list);
    virtual ~GCTAEventList(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCTAEventList* clone(void) const;
    virtual int            size(void) const { return m_events.size(); }
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename, bool clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const { return m_events.size(); }
    virtual void           roi(const GRoi& roi);
    virtual const GCTARoi& roi(void) const { return m_roi; }

    // Implement other methods
    void                   append(const GCTAEventAtom& event);
    void                   reserve(const int& number);
};


/***********************************************************************//**
 * @brief GCTAEventList class extension
 ***************************************************************************/
%extend GCTAEventList {
    GCTAEventList copy() {
        return (*self);
    }
    GCTAEventAtom* __getitem__(int index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GCTAEventAtom& val) {
        if (index>=0 && index < self->size())
            *((*self)[index]) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
};


/***********************************************************************//**
 * @brief GCTAEventList type casts
 ***************************************************************************/
%inline %{
    GCTAEventList* cast_GCTAEventList(GEvents* events) {
        GCTAEventList* list = dynamic_cast<GCTAEventList*>(events);
        if (list == NULL)
            throw GException::bad_type("cast_GCTAEventList(GEvents*)",
                                       "GEvents not of type GCTAEventList");            
        return list;
    }
%}
