/***************************************************************************
 *      GModelSpatialGauss.i  -  Spatial Gaussian source model class       *
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
 * @file GModelSpatialGauss.i
 * @brief Gaussian spatial model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialGauss.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialGauss
 *
 * @brief Gaussian spatial model class
 ***************************************************************************/
class GModelSpatialGauss : public GModelSpatial {
public:
    // Constructors and destructors
    GModelSpatialGauss(void);
    explicit GModelSpatialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelSpatialGauss(const GXmlElement& xml);
    GModelSpatialGauss(const GModelSpatialGauss& model);
    virtual ~GModelSpatialGauss(void);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialGauss* clone(void) const;
    virtual std::string         type(void) const;
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;

    // Other methods
    double  ra(void) const;
    double  dec(void) const;
    double  sigma(void) const;
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    sigma(const double& sigma);
};


/***********************************************************************//**
 * @brief GModelSpatialGauss class extension
 ***************************************************************************/
%extend GModelSpatialGauss {
    GModelSpatialGauss copy() {
        return (*self);
    }
};
