/***************************************************************************
 *     GModelRadialRegistry.i  -  Radial spatial model registry class      *
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
 * @file GModelRadialRegistry.i
 * @brief Radial spatial model registry class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelRadialRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelRadialRegistry
 *
 * @brief Radial spatial model registry class
 ***************************************************************************/
class GModelRadialRegistry {
public:
    // Constructors and destructors
    GModelRadialRegistry(void);
    explicit GModelRadialRegistry(const GModelRadial* model);
    GModelRadialRegistry(const GModelRadialRegistry& registry);
    virtual ~GModelRadialRegistry(void);

    // Methods
    int           size(void) const;
    GModelRadial* alloc(const std::string& type) const;
    std::string   name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelRadialRegistry class extension
 ***************************************************************************/
%extend GModelRadialRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
