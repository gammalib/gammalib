/***************************************************************************
 *           GInstDir.i  -  Instrument direction class python I/F          *
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
 * @file GInstDir.i
 * @brief GInstDir class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GInstDir.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GInstDir
 *
 * @brief Abstract interface for the instrument direction of an event.
 *
 * The instrument direction of an event is the equivalent of the sky
 * direction (implemented by GSkyDir) but in the instrument data space.
 * The instrument direction may be any kind of position or direction
 * information encoded in the data space, such as incident event
 * reconstructions for imaging devices or detector numbers etc. for
 * non-imaging devices.
 ***************************************************************************/
class GInstDir {
public:
    // Constructors and destructors
    GInstDir(void);
    GInstDir(const GInstDir& dir);
    virtual ~GInstDir(void);

    // Pure virtual methods
    virtual void      clear(void) = 0;
    virtual GInstDir* clone(void) const = 0;
};


/***********************************************************************//**
 * @brief GInstDir class extension
 ***************************************************************************/
%extend GInstDir {
    char *__str__() {
        return tochar(self->print());
    }
};
