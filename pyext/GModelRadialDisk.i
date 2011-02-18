/***************************************************************************
 *        GModelRadialDisk.i  -  Radial disk source model class            *
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
 * @file GModelRadialDisk.i
 * @brief Radial disk model class Python interface definition
 * @author C. Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelRadialDisk.hpp"
#include "GTools.hpp"
%}

/**************************************************************************
 * @class GModelRadialDisk
 *
 * @brief Disk source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a disk source, i.e. constant surface brightness within some
 * radius and no emission outside.
 ***************************************************************************/
class GModelRadialDisk : public GModelRadial {

public:
    // Constructors and destructors
    GModelRadialDisk(void);
    explicit GModelRadialDisk(const GSkyDir& dir, const double& radius);
    explicit GModelRadialDisk(const GXmlElement& xml);
    GModelRadialDisk(const GModelRadialDisk& model);
    virtual ~GModelRadialDisk(void);

    // Implemented pure virtual methods
    virtual void              clear(void);
    virtual GModelRadialDisk* clone(void) const;
    virtual std::string       type(void) const;
    virtual double            eval(const double& theta) const;
    virtual double            eval_gradients(const double& theta) const;
    virtual GSkyDir           mc(GRan& ran) const;
    virtual double            theta_max(void) const;
    virtual void              read(const GXmlElement& xml);
    virtual void              write(GXmlElement& xml) const;

    // Other methods
    double radius(void) const;
    void   radius(const double& radius);
};


/***********************************************************************//**
 * @brief GModelRadialDisk class extension
 ***************************************************************************/
%extend GModelRadialDisk {
    GModelRadialDisk copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GModelRadialDisk type casts
 ***************************************************************************/
%inline %{
    GModelRadialDisk* cast_GModelRadialDisk(GModelSpatial* model) {
        return dynamic_cast<GModelRadialDisk*>(model);
    }
%};
