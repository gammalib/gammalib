/***************************************************************************
 *          GModelSpatial.i  -  Spatial model class SWIG interface         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatial.i
 * @brief GModelSpatial class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatial.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatial
 *
 * @brief Abstract SWIG interface definition for the spatial model class.
 ***************************************************************************/
class GModelSpatial {
public:
    // Constructors and destructors
    explicit GModelSpatial(void);
    GModelSpatial(const GModelSpatial& model);
    virtual ~GModelSpatial();

    // Virtual methods
    virtual int        npars(void) const = 0;
    virtual GModelPar* par(int index) const = 0;
    virtual double     eval(const GSkyDir& srcDir) = 0;
    virtual double     eval_gradients(const GSkyDir& srcDir) = 0;
    virtual bool       isptsource(void) const { return false; }
};
