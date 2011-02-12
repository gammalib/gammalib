/***************************************************************************
 *      GModelSpatialShell.i  -  Spatial shell source model class          *
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
 * @file GModelSpatialShell.i
 * @brief Spatial shell model class Python interface definition
 * @author C. Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialShell.hpp"
#include "GTools.hpp"
%}

/**************************************************************************
 * @class GModelSpatialShell
 *
 * @brief Shell source model class
 ***************************************************************************/
class GModelSpatialShell : public GModelSpatial {
	public:
    // Constructors and destructors
    GModelSpatialShell(void);
    explicit GModelSpatialShell(const GSkyDir& dir,
								const double& theta_in, const double& theta_out);
    explicit GModelSpatialShell(const GXmlElement& xml);
    GModelSpatialShell(const GModelSpatialShell& model);
    virtual ~GModelSpatialShell(void);
	
    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialShell* clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    double  ra(void);
    double  dec(void);
    double  radius(void) const;
    double  width(void) const;
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    radius(const double& radius);
    void    width(const double& width);
};


/***********************************************************************//**
 * @brief GModelSpatialShell class extension
 ***************************************************************************/
%extend GModelSpatialShell {
    GModelSpatialShell copy() {
        return (*self);
    }
};
