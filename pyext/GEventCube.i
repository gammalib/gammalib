/***************************************************************************
 *     GEventCube.i  -  Abstract event cube container class python I/F     *
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
 * @file GEventCube.i
 * @brief GEventCube class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEventCube.hpp"
%}


/***********************************************************************//**
 * @class GEventCube
 *
 * @brief GEventCube container class interface defintion.
 ***************************************************************************/
class GEventCube : public GEvents {
public:
    // Constructors and destructors
    GEventCube(void);
    GEventCube(const GEventCube& cube);
    virtual ~GEventCube(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventCube* clone(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual GEventBin*  pointer(int index) = 0;
    virtual int         number(void) const = 0;

    // Implemented pure virtual methods
    int  size(void) const { return m_elements; }
    int  dim(void) const { return m_dim; }
    int  naxis(int axis) const;
    bool islist(void) const { return false; }
    bool iscube(void) const { return true; }
};


/***********************************************************************//**
 * @brief GEventCube class extension
 *
 * The __getitem__ method makes the event cube iteratable.
 ***************************************************************************/
%extend GEventCube {
    GEventBin* __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return self->pointer(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
};


/***********************************************************************//**
 * @brief GEventCube type casts
 ***************************************************************************/
%inline %{
    GEventCube* cast_GEventCube(GEvents* events) {
        if (!events->iscube())
            throw GException::fits_invalid_type("cast_GEventCube(GEvents*)",
                                                "GEvents is not an event cube.");
        return dynamic_cast<GEventCube*>(events);
    }
%};
