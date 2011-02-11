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
 * @brief Radial Gaussian model class Python interface definition
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
 * @brief Radial Gaussian CTA model class
 ***************************************************************************/
class GCTAModelRadialGauss : public GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadialGauss(void);
    explicit GCTAModelRadialGauss(const double& sigma);
    explicit GCTAModelRadialGauss(const GXmlElement& xml);
    GCTAModelRadialGauss(const GCTAModelRadialGauss& model);
    virtual ~GCTAModelRadialGauss(void);

    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GCTAModelRadialGauss* clone(void) const;
    virtual std::string           type(void) const;
    virtual double                eval(const double& offset) const;
    virtual double                eval_gradients(const double& offset) const;
    virtual GCTAInstDir           mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                omega(void) const;
    virtual void                  read(const GXmlElement& xml);
    virtual void                  write(GXmlElement& xml) const;

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
};
