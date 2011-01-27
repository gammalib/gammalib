/***************************************************************************
 * GCTAModelRadialGauss.i  -  Radial Gaussian model class python interface *
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
 * @file GCTAModelRadialGauss.i
 * @brief GCTAModelRadialGauss class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialGauss.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialGauss
 *
 * @brief Radial Gaussian CTA model interface definition
 ***************************************************************************/
class GCTAModelRadialGauss  : public GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadialGauss(void);
    explicit GCTAModelRadialGauss(const double& sigma);
    explicit GCTAModelRadialGauss(const GXmlElement& xml);
    GCTAModelRadialGauss(const GCTAModelRadialGauss& model);
    virtual ~GCTAModelRadialGauss(void);

    // Implemented pure virtual methods
    void                  clear(void);
    GCTAModelRadialGauss* clone(void) const;
    int                   size(void) const;
    std::string           type(void) const;
    double                eval(const double& offset);
    double                eval_gradients(const double& offset);
    GCTAInstDir           mc(const GCTAInstDir& dir, GRan& ran) const;
    double                omega(void) const;
    void                  read(const GXmlElement& xml);
    void                  write(GXmlElement& xml) const;

    // Other methods
    double sigma(void) const;
    void   sigma(const double& sigma);
};


/***********************************************************************//**
 * @brief GCTAModelRadialGauss class extension
 ***************************************************************************/
%extend GCTAModelRadialGauss {
    char *__str__() {
        return tochar(self->print());
    }
    /*
    GModelPar __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GModelPar& val) {
        if (index>=0 && index < self->size())
            (*self)(index) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    */
};
