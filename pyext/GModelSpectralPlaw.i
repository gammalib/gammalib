/***************************************************************************
 * GModelSpectralPlaw.i  -  Spectral power law model class SWIG interface  *
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
 * @file GModelSpectralPlaw.i
 * @brief GModelSpectralPlaw class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralPlaw.hpp"
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
    GModelSpectralPlaw(const GModelSpectralPlaw& model);
    virtual ~GModelSpectralPlaw(void);

    // Methods
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
    GModelPar* par_norm(void) { return &m_norm; }
    GModelPar* par_index(void) { return &m_index; }
    GModelPar* par_pivot(void) { return &m_pivot; }
    double     eval(const GEnergy& srcEng);
    double     eval_gradients(const GEnergy& srcEng);
    void       autoscale(void);
    double     norm(void) const { return m_norm.real_value(); }
    double     index(void) const { return m_index.real_value(); }
    double     pivot(void) const { return m_pivot.real_value(); }
};


/***********************************************************************//**
 * @brief GModelSpectralPlaw class extension
 ***************************************************************************/
%extend GModelSpectralPlaw {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    GModelSpectralPlaw copy() {
        return (*self);
    }
};
