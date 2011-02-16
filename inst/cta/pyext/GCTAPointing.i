/***************************************************************************
 *                   GCTAPointing.i  -  CTA pointing class                 *
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
 * @file GCTAPointing.i
 * @brief CTA pointing class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAPointing.hpp"
%}


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief CTA pointing class
 ***************************************************************************/
class GCTAPointing : public GPointing {
public:
    // Constructors and destructors
    GCTAPointing(void);
    GCTAPointing(const GCTAPointing& pnt);
    ~GCTAPointing(void);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GCTAPointing*  clone(void) const;
    virtual const GSkyDir& dir(void) const;

    // Other methods
    void dir(const GSkyDir& dir);
};


/***********************************************************************//**
 * @brief GCTAPointing class extension
 ***************************************************************************/
%extend GCTAPointing {
    GCTAPointing copy() {
        return (*self);
    }
};
