/***************************************************************************
 *              GPointing.i  -  Pointing class python bindings             *
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
 * @file GPointing.i
 * @brief GPointing class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPointing.hpp"
%}


/***********************************************************************//**
 * @class GPointing
 *
 * @brief Pointing class
 ***************************************************************************/
class GPointing {

public:
    // Constructors and destructors
    GPointing(void);
    GPointing(const GPointing& pnt);
    virtual ~GPointing(void);

    // Pure virtual methods
    virtual GPointing* clone(void) const = 0;
};
