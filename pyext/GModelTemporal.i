/***************************************************************************
 *         GModelTemporal.i  -  Temporal model class SWIG interface        *
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
 * @file GModelTemporal.i
 * @brief GModelTemporal class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporal.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporal
 *
 * @brief Abstract SWIG interface definition for the temporal model class.
 ***************************************************************************/
class GModelTemporal {
public:
    // Constructors and destructors
    GModelTemporal(void);
    GModelTemporal(const GModelTemporal& model);
    virtual ~GModelTemporal(void);

    // Virtual methods
    virtual void            clear(void) = 0;
    virtual GModelTemporal* clone(void) const = 0;
    virtual int             size(void) const = 0;
    virtual std::string     type(void) const = 0;
    virtual GModelPar*      par(int index) const = 0;
    virtual double          eval(const GTime& srcTime) = 0;
    virtual double          eval_gradients(const GTime& srcTime) = 0;
};
