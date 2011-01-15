/***************************************************************************
 * GModelSpectralPlaw.i  -  Spectral power law model class SWIG interface  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralPlaw.i
 * @brief GModelSpectralPlaw class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlaw.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralPlaw
 *
 * @brief Powerlaw SWIG interface definition.
 ***************************************************************************/
class GModelSpectralPlaw  : public GModelSpectral {
public:
    // Constructors and destructors
    GModelSpectralPlaw(void);
    explicit GModelSpectralPlaw(const double& norm, const double& index);
    explicit GModelSpectralPlaw(const GXmlElement& xml);
    GModelSpectralPlaw(const GModelSpectralPlaw& model);
    virtual ~GModelSpectralPlaw(void);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpectralPlaw* clone(void) const;
    int                 size(void) const { return m_npars; }
    std::string         type(void) const { return "PowerLaw"; }
    double              eval(const GEnergy& srcEng);
    double              eval_gradients(const GEnergy& srcEng);
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;

    // Other methods
    void   autoscale(void);
    double norm(void) const { return m_norm.real_value(); }
    double index(void) const { return m_index.real_value(); }
    double pivot(void) const { return m_pivot.real_value(); }
};


/***********************************************************************//**
 * @brief GModelSpectralPlaw class extension
 ***************************************************************************/
%extend GModelSpectralPlaw {
    char *__str__() {
        return tochar(self->print());
    }
    GModelSpectralPlaw copy() {
        return (*self);
    }
};
