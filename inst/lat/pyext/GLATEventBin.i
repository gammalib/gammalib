/***************************************************************************
 *                  GLATEventBin.i  -  LAT event bin class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventBin.i
 * @brief GLATEventBin class python bindings.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventBin.hpp"
%}


/***********************************************************************//**
 * @class GLATEventBin
 *
 * @brief GLATEventBin class interface defintion.
 ***************************************************************************/
class GLATEventBin : public GEventBin {

    // Friend classes
    friend class GLATEventCube;

public:
    // Constructors and destructors
    GLATEventBin(void);
    GLATEventBin(const GLATEventBin& bin);
    virtual ~GLATEventBin(void);

    // Implemented pure virtual base class methods
    void               clear(void);
    GLATEventBin*      clone(void) const;
    double             size(void) const;
    const GLATInstDir& dir(void) const { return *m_dir; }
    const GEnergy&     energy(void) const { return *m_energy; }
    const GTime&       time(void) const { return *m_time; }
    double             counts(void) const { return *m_counts; }
    double             error(void) const;

    // Other methods
    const double&  omega(void) const { return *m_omega; }
    const GEnergy& ewidth(void) const { return *m_ewidth; }
    const double&  ontime(void) const { return *m_ontime; }
    const int&     index(void) const { return m_index; }
    const int&     ipix(void) const { return m_ipix; }
    const int&     ieng(void) const { return m_ieng; }
    GLATEventCube* cube(void) const { return m_cube; }
};


/***********************************************************************//**
 * @brief GLATEventBin class extension
 ***************************************************************************/
%extend GLATEventBin {
    GLATEventBin(const GEvent& event) {
        if (!event.isbin())
            throw GException::bad_type("GLATEventBin(GEvent&)",
                                       "GEvent not an event bin");            
        GLATEventBin* bin = new GLATEventBin((GLATEventBin&)event);
        return bin;
    }
    GLATEventBin(const GEventBin& event) {
        GLATEventBin* bin = new GLATEventBin((GLATEventBin&)event);
        return bin;
    }
    GLATEventBin copy() {
        return (*self);
    }
};
