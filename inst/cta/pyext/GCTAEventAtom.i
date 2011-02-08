/***************************************************************************
 *                 GCTAEventAtom.i  -  CTA event atom class                *
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
 * @file GCTAEventAtom.i
 * @brief CTA event bin class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventAtom.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventAtom
 *
 * @brief CTA event atom class Python interface
 ***************************************************************************/
class GCTAEventAtom : public GEventAtom {
public:
    // Constructors and destructors
    GCTAEventAtom(void);
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom(void);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCTAEventAtom*     clone(void) const;
    const GCTAInstDir& dir(void) const { return m_dir; }
    const GEnergy&     energy(void) const { return m_energy; }
    const GTime&       time(void) const { return m_time; }
    void               dir(const GCTAInstDir& dir) { m_dir=dir; }
    void               energy(const GEnergy& energy) { m_energy=energy; }
    void               time(const GTime& time) { m_time=time; }
};


/***********************************************************************//**
 * @brief GCTAEventAtom class extension
 ***************************************************************************/
%extend GCTAEventAtom {
    GCTAEventAtom copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GCTAEventAtom type casts
 ***************************************************************************/
%inline %{
    GCTAEventAtom* cast_GCTAEventAtom(GEvent* event) {
        GCTAEventAtom* atom = dynamic_cast<GCTAEventAtom*>(event);
        if (atom == NULL)
            throw GException::bad_type("cast_GCTAEventAtom(GEvent*)",
                                       "GEvent not of type GCTAEventAtom");            
        return atom;
    }
%}
