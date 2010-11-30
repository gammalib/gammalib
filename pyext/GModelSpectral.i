/***************************************************************************
 *         GModelSpectral.i  -  Spectral model class SWIG interface        *
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
 * @file GModelSpectral.i
 * @brief GModelSpectral class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectral.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectral
 *
 * @brief Abstract SWIG interface definition for the spectral model class.
 ***************************************************************************/
class GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectral(void);
    GModelSpectral(const GModelSpectral& model);
    virtual ~GModelSpectral();

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GModelSpectral* clone(void) const = 0;
    virtual int             size(void) const = 0;
    virtual std::string     name(void) const = 0;
    virtual GModelPar*      par(int index) const = 0;
    virtual double          eval(const GEnergy& srcEng) = 0;
    virtual double          eval_gradients(const GEnergy& srcEng) = 0;
    virtual void            read(const GXmlElement& xml) = 0;
    virtual void            write(GXmlElement& xml) const = 0;
};
