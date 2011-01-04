/***************************************************************************
 *                 GLATEventCube.i  -  LAT event cube class                *
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
 * @file GLATEventCube.i
 * @brief GLATEventCube class python bindings.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventCube.hpp"
%}


/***********************************************************************//**
 * @class GLATEventCube
 *
 * @brief GLATEventCube class interface defintion.
 ***************************************************************************/
class GLATEventCube : public GEventCube {

public:
    // Constructors and destructors
    GLATEventCube(void);
    GLATEventCube(const GLATEventCube& cube);
    virtual ~GLATEventCube(void);

    // Implemented pure virtual base class methods
    void           clear(void);
    GLATEventCube* clone(void) const;
    void           load(const std::string& filename);
    GLATEventBin*  pointer(int index);
    int            number(void) const;

    // Other methods
    GEbounds&   ebds(void) { return m_ebds; }
    GGti&       gti(void) { return m_gti; }
    GTime&      time(void) { return m_time; }
    GSkymap&    map(void) { return m_map; }
    double      ontime(void) { return m_ontime; }
    int         nx(void) const { return m_map.nx(); }
    int         ny(void) const { return m_map.ny(); }
    int         npix(void) const { return m_map.npix(); }
    int         ebins(void) const { return m_map.nmaps(); }
    int         ndiffrsp(void) const { return m_srcmap.size(); }
    std::string diffname(const int& index) const;
    GSkymap*    diffrsp(const int& index) const;
    GNodeArray* enodes(void) { return &m_enodes; }
    double      maxrad(const GSkyDir& dir) const;
};


/***********************************************************************//**
 * @brief GLATEventCube class extension
 *
 * The GLATEventCube methods perform type conversion.
 * The __getitem__ method makes the event cube iteratable.
 ***************************************************************************/
%extend GLATEventCube {
    GLATEventCube(const GEvents& events) {
        if (!events.iscube())
            throw GException::bad_type("GLATEventCube(GEvents&)",
                                       "GEvents not an event cube");            
        GLATEventCube* cube = (GLATEventCube*)&events;
        return cube;
    }
    GLATEventCube(const GEventCube& events) {
        GLATEventCube* cube = (GLATEventCube*)&events;
        return cube;
    }
    GLATEventBin* __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return self->pointer(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    GLATEventCube copy() {
        return (*self);
    }
};
