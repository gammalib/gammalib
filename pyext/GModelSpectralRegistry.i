/***************************************************************************
 *  GModelSpectralRegistry.i  -  Spectral model registry class python I/F  *
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
 * @file GModelSpectralRegistry.i
 * @brief GModelSpectralRegistry class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpectralRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpectralRegistry
 *
 * @brief Interface definition for the spectral model registry class.
 ***************************************************************************/
class GModelSpectralRegistry {
public:
    // Constructors and destructors
    GModelSpectralRegistry(void);
    GModelSpectralRegistry(const GModelSpectral* model);
    GModelSpectralRegistry(const GModelSpectralRegistry& registry);
    virtual ~GModelSpectralRegistry(void);

    // Methods
    int             size(void) const { return m_number; }
    GModelSpectral* alloc(const std::string& type) const;
    std::string     name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelSpectralRegistry class extension
 ***************************************************************************/
%extend GModelSpectralRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
