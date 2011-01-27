/***************************************************************************
 * GCTAModelRadialRegistry.i  -  CTA Radial model registry class python I/F*
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
 * @file GCTAModelRadialRegistry.i
 * @brief GCTAModelRadialRegistry class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialRegistry
 *
 * @brief Interface definition for the spatial model registry class.
 ***************************************************************************/
class GCTAModelRadialRegistry {
public:
    // Constructors and destructors
    GCTAModelRadialRegistry(void);
    GCTAModelRadialRegistry(const GCTAModelRadial* model);
    GCTAModelRadialRegistry(const GCTAModelRadialRegistry& registry);
    virtual ~GCTAModelRadialRegistry(void);

    // Methods
    int              size(void) const;
    GCTAModelRadial* alloc(const std::string& type) const;
    std::string      name(const int& index) const;
};


/***********************************************************************//**
 * @brief GCTAModelRadialRegistry class extension
 ***************************************************************************/
%extend GCTAModelRadialRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
