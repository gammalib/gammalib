/***************************************************************************
 *               GCTARoi.i  -  CTA region of interest class                *
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
 * @file GCTARoi.i
 * @brief CTA region of interest class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTARoi.hpp"
%}


/***********************************************************************//**
 * @class GCTARoi
 *
 * @brief CTA region of interest class
 ***************************************************************************/
class GCTARoi : public GRoi {

public:
    // Constructors and destructors
    GCTARoi(void);
    GCTARoi(const GCTARoi& roi);
    virtual ~GCTARoi(void);

    // Implemented pure virtual base class methods
    void        clear(void);
    GCTARoi*    clone(void) const;

    // Other methods
    GCTAInstDir centre(void) const;
    double      radius(void) const;
    void        centre(const GCTAInstDir& centre);
    void        radius(const double& radius);
};


/***********************************************************************//**
 * @brief GCTARoi class extension
 ***************************************************************************/
%extend GCTARoi {
    GCTARoi copy() {
        return (*self);
    }
};
