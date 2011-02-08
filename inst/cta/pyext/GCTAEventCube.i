/***************************************************************************
 *             GCTAEventCube.i  -  CTA event bin container class           *
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
 * @file GCTAEventCube.i
 * @brief CTA event bin container class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAEventCube.hpp"
%}


/***********************************************************************//**
 * @class GCTAEventCube
 *
 * @brief CTA event bin container class Python interface
 ***************************************************************************/
class GCTAEventCube : public GEventCube {
public:
    // Constructors and destructors
    GCTAEventCube(void);
    explicit GCTAEventCube(const GSkymap& map, const GEbounds& ebds, const GGti& gti);
    GCTAEventCube(const GCTAEventCube& cube);
    virtual ~GCTAEventCube(void);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCTAEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(int axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename, bool clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;

    // Other methods
    void                   map(const GSkymap& map) { m_map=map; }
    const GSkymap&         map(void) const { return m_map; }
    int                    nx(void) const { return m_map.nx(); }
    int                    ny(void) const { return m_map.ny(); }
    int                    npix(void) const { return m_map.npix(); }
    int                    ebins(void) const { return m_map.nmaps(); }
};


/***********************************************************************//**
 * @brief GCTAEventCube class extension
 ***************************************************************************/
%extend GCTAEventCube {
    GCTAEventCube copy() {
        return (*self);
    }
    GCTAEventBin* __getitem__(int index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GCTAEventBin& val) {
        if (index>=0 && index < self->size())
            *((*self)[index]) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
};


/***********************************************************************//**
 * @brief GCTAEventCube type casts
 ***************************************************************************/
%inline %{
    GCTAEventCube* cast_GCTAEventCube(GEvents* events) {
        GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>(events);
        if (cube == NULL)
            throw GException::bad_type("cast_GCTAEventCube(GEvents*)",
                                       "GEvents not of type GCTAEventCube");            
        return cube;
    }
%}
