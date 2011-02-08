/***************************************************************************
 *                  GLATEventBin.i  -  LAT event bin class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @brief LAT event bin class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATEventBin.hpp"
%}


/***********************************************************************//**
 * @class GLATEventBin
 *
 * @brief LAT event bin class Python interface
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
    const GLATInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    double             counts(void) const;
    double             error(void) const;

    // Other methods
    const double&  omega(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
    const int&     index(void) const;
    const int&     ipix(void) const;
    const int&     ieng(void) const;
    GLATEventCube* cube(void) const;
};


/***********************************************************************//**
 * @brief GLATEventBin class extension
 ***************************************************************************/
%extend GLATEventBin {
    GLATEventBin copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GLATEventBin type casts
 ***************************************************************************/
%inline %{
    GLATEventBin* cast_GLATEventBin(GEvent* event) {
        GLATEventBin* bin = dynamic_cast<GLATEventBin*>(event);
        if (bin == NULL)
            throw GException::bad_type("cast_GLATEventBin(GEvent*)",
                                       "GEvent not of type GLATEventBin");            
        return bin;
    }
%}
