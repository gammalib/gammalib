/***************************************************************************
 *                 GPointing.i  -  Abstract pointing class                 *
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
 * @file GPointing.i
 * @brief Abstract pointing class Python interface definition.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPointing.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPointing
 *
 * @brief Abstract Pointing class
 ***************************************************************************/
class GPointing {

public:
    // Constructors and destructors
    GPointing(void);
    GPointing(const GPointing& pnt);
    virtual ~GPointing(void);

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GPointing*     clone(void) const = 0;
    virtual const GSkyDir& dir(void) const = 0;
};


/***********************************************************************//**
 * @brief GPointing class extension
 ***************************************************************************/
%extend GPointing {
    char *__str__() {
        return tochar(self->print());
    }
};
