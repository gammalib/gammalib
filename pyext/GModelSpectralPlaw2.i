/***************************************************************************
 *   GModelSpectralPlaw2.i  -  Spectral power law model class python I/F   *
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
 * @file GModelSpectralPlaw2.i
 * @brief GModelSpectralPlaw2 class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlaw2.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlaw2
 *
 * @brief Powerlaw python interface definition
 ***************************************************************************/
class GModelSpectralPlaw2  : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlaw2(void);
    explicit GModelSpectralPlaw2(const double& integral, const double& index);
    explicit GModelSpectralPlaw2(const GXmlElement& xml);
    GModelSpectralPlaw2(const GModelSpectralPlaw2& model);
    virtual ~GModelSpectralPlaw2(void);

    // Implemented pure virtual methods
    void                 clear(void);
    GModelSpectralPlaw2* clone(void) const;
    int                  size(void) const;
    std::string          type(void) const;
    double               eval(const GEnergy& srcEng);
    double               eval_gradients(const GEnergy& srcEng);
    double               flux(const GEnergy& emin, const GEnergy& emax) const;
    GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    void                 read(const GXmlElement& xml);
    void                 write(GXmlElement& xml) const;

    // Other methods
    double integral(void) const;
    double index(void) const;
    double emin(void) const;
    double emax(void) const;
};


/***********************************************************************//**
 * @brief GModelSpectralPlaw2 class extension
 ***************************************************************************/
%extend GModelSpectralPlaw2 {
    GModelSpectralPlaw2 copy() {
        return (*self);
    }
};
