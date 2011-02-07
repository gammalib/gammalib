/***************************************************************************
 *                 GLATEventCube.i  -  LAT event cube class                *
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
    void              clear(void);
    GLATEventCube*    clone(void) const;
    int               size(void) const;
    int               dim(void) const;
    int               naxis(int axis) const;
    void              load(const std::string& filename);
    GLATEventBin*     pointer(int index);
    int               number(void) const;

    // Other methods
    void              ebds(const GEbounds& ebds);
    void              gti(const GGti& gti);
    void              time(const GTime& time);
    void              map(const GSkymap& map);
    void              enodes(const GNodeArray& enodes);
    void              ontime(const double& ontime);
    const GEbounds&   ebds(void) const;
    const GGti&       gti(void) const;
    const GTime&      time(void) const;
    const GSkymap&    map(void) const;
    const GNodeArray& enodes(void);
    const double&     ontime(void) const;
    int               nx(void) const;
    int               ny(void) const;
    int               npix(void) const;
    int               ebins(void) const;
    int               ndiffrsp(void) const;
    std::string       diffname(const int& index) const;
    GSkymap*          diffrsp(const int& index) const;
    double            maxrad(const GSkyDir& dir) const;
};


/***********************************************************************//**
 * @brief GLATEventCube class extension
 *
 * The __getitem__ method makes the event cube iteratable.
 ***************************************************************************/
%extend GLATEventCube {
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


/***********************************************************************//**
 * @brief GLATEventCube type casts
 ***************************************************************************/
%inline %{
    GLATEventCube* cast_GLATEventCube(GEvents* events) {
        if (!events->iscube())
            throw GException::bad_type("cast_GLATEventCube(GEvents*)",
                                       "GEvents not an event cube");            
        return dynamic_cast<GLATEventCube*>(events);
    }
    GLATEventCube* cast_GLATEventCube(GEventCube* events) {
        return dynamic_cast<GLATEventCube*>(events);
    }
%}
