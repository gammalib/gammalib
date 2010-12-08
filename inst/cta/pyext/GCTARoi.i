/***************************************************************************
 *               GCTARoi.i  -  CTA region of interest class                *
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
 * @file GCTARoi.i
 * @brief GCTARoi class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTARoi.hpp"
%}


/***********************************************************************//**
 * @class GCTARoi
 *
 * @brief Python bindings for the CTA region of interest class
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
    GCTAInstDir centre(void) const { return m_centre; }
    double      radius(void) const { return m_radius; }
    void        centre(const GCTAInstDir& centre) { m_centre=centre; return; }
    void        radius(const double& radius) { m_radius=radius; return; }
};


/***********************************************************************//**
 * @brief GCTARoi class extension
 ***************************************************************************/
%extend GCTARoi {
    GCTARoi copy() {
        return (*self);
    }
};
