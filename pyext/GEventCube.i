/***************************************************************************
 *           GEventCube.i  -  Abstract event bin container class           *
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
 * @brief Abstract event bin container class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEventCube.hpp"
%}


/***********************************************************************//**
 * @class GEventCube
 *
 * @brief Abstract event bin container class interface
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
    virtual int         size(void) const = 0;
    virtual int         dim(void) const = 0;
    virtual int         naxis(int axis) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual void        save(const std::string& filename, bool clobber = false) const = 0;
    virtual void        read(const GFits& file) = 0;
    virtual void        write(GFits& file) const = 0;
    virtual int         number(void) const = 0;
};


/***********************************************************************//**
 * @brief GEventCube class extension
 ***************************************************************************/
%extend GEventCube {
};


/***********************************************************************//**
 * @brief GEventCube type casts
 ***************************************************************************/
%inline %{
    GEventCube* cast_GEventCube(GEvents* events) {
        GEventCube* cube = dynamic_cast<GEventCube*>(events);
        if (cube == NULL)
            throw GException::fits_invalid_type("cast_GEventCube(GEvents*)",
                                                "GEvents is not of type GEventCube.");
        return cube;
    }
%};
