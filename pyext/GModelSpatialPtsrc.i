/***************************************************************************
 *       GModelSpatialPtsrc.i  -  Spatial point source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialPtsrc.i
 * @brief Point source spatial model class Python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialPtsrc.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialPtsrc
 *
 * @brief Point source spatial model
 ***************************************************************************/
class GModelSpatialPtsrc  : public GModelSpatial {
public:
    // Constructors and destructors
    explicit GModelSpatialPtsrc(void);
    explicit GModelSpatialPtsrc(const GSkyDir& dir);
    explicit GModelSpatialPtsrc(const GXmlElement& xml);
    GModelSpatialPtsrc(const GModelSpatialPtsrc& model);
    virtual ~GModelSpatialPtsrc(void);

    // Implemented virtual methods
    virtual void                clear(void);
    virtual GModelSpatialPtsrc* clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    double  ra(void) const;
    double  dec(void) const;
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
};


/***********************************************************************//**
 * @brief GModelSpatialPtsrc class extension
 ***************************************************************************/
%extend GModelSpatialPtsrc {
    GModelSpatialPtsrc copy() {
        return (*self);
    }
};
