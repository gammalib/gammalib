/***************************************************************************
 *                  GCTAEventBin.i  -  CTA event bin class                 *
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
 * @file GCTAEventBin.i
 * @brief CTA event bin class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventBin.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventBin
 *
 * @brief CTA event bin class Python interface
 ***************************************************************************/
class GCTAEventBin : public GEventBin {

    // Friend classes
    friend class GCTAEventCube;

public:
    // Constructors and destructors
    GCTAEventBin(void);
    GCTAEventBin(const GCTAEventBin& bin);
    virtual ~GCTAEventBin(void);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCTAEventBin*      clone(void) const;
    double             size(void) const;
    const GCTAInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    double             counts(void) const;
    double             error(void) const;

    // Other methods
    const double&  omega(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
};


/***********************************************************************//**
 * @brief GCTAEventBin class extension
 ***************************************************************************/
%extend GCTAEventBin {
    GCTAEventBin copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GCTAEventBin type casts
 ***************************************************************************/
%inline %{
    GCTAEventBin* cast_GCTAEventBin(GEvent* event) {
        GCTAEventBin* bin = dynamic_cast<GCTAEventBin*>(event);
        if (bin == NULL)
            throw GException::bad_type("cast_GCTAEventBin(GEvent*)",
                                       "GEvent not of type GCTAEventBin");            
        return bin;
    }
%}
