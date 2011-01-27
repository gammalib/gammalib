/***************************************************************************
 * GModelSpatialGauss.i  -  Spatial Gaussian source model class python I/F *
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
 * @brief GModelSpatialGauss class python interface.
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
 * @brief Gaussian source model interface definition.
 ***************************************************************************/
class GModelSpatialGauss  : public GModelSpatial {
public:
    // Constructors and destructors
    GModelSpatialGauss(void);
    explicit GModelSpatialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelSpatialGauss(const GXmlElement& xml);
    GModelSpatialGauss(const GModelSpatialGauss& model);
    virtual ~GModelSpatialGauss(void);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpatialGauss* clone(void) const;
    int                 size(void) const;
    std::string         type(void) const;
    double              eval(const GSkyDir& srcDir);
    double              eval_gradients(const GSkyDir& srcDir);
    GSkyDir             mc(GRan& ran) const;
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    bool                isptsource(void) const;

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
