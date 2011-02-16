/***************************************************************************
 *        GModelSpatialDisk.i  -  Spatial disk source model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Christoph Deil                                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialDisk.i
 * @brief Disk spatial model class Python interface definition
 * @author C. Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialDisk.hpp"
#include "GTools.hpp"
%}

/**************************************************************************
 * @class GModelSpatialDisk
 *
 * @brief Disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a disk source, i.e. constant surface brightness within some
 * radius and no emission outside.
 ***************************************************************************/
class GModelSpatialDisk : public GModelSpatial {
	
public:
    // Constructors and destructors
    GModelSpatialDisk(void);
    explicit GModelSpatialDisk(const GSkyDir& dir,
							   const double& radius);
    explicit GModelSpatialDisk(const GXmlElement& xml);
    GModelSpatialDisk(const GModelSpatialDisk& model);
    virtual ~GModelSpatialDisk(void);
	
    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialDisk*  clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
	
    // Other methods
    double  ra(void) const;
    double  dec(void) const;
    double  radius(void) const;
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    radius(const double& radius);
};


/***********************************************************************//**
 * @brief GModelSpatialDisk class extension
 ***************************************************************************/
%extend GModelSpatialDisk {
    GModelSpatialDisk copy() {
        return (*self);
    }
};
