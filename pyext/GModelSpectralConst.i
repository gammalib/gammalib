/***************************************************************************
 *         GModelSpectralConst.i  -  Spectral constant model class         *
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
 * @file GModelSpectralConst.i
 * @brief Spectral contant model class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralConst.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralConst
 *
 * @brief Spectral constant model class
 ***************************************************************************/
class GModelSpectralConst : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralConst(void);
    explicit GModelSpectralConst(const GXmlElement& xml);
    GModelSpectralConst(const GModelSpectralConst& model);
    virtual ~GModelSpectralConst(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GModelSpectralConst* clone(void) const;
    virtual std::string          type(void) const;
    virtual double               eval(const GEnergy& srcEng) const;
    virtual double               eval_gradients(const GEnergy& srcEng) const;
    virtual double               flux(const GEnergy& emin, const GEnergy& emax) const;
    virtual GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;

    // Other methods
    double norm(void) const;
};


/***********************************************************************//**
 * @brief GModelSpectralConst class extension
 ***************************************************************************/
%extend GModelSpectralConst {
    GModelSpectralConst copy() {
        return (*self);
    }
};
