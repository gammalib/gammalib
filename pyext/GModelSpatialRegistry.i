/***************************************************************************
 *   GModelSpatialRegistry.i  -  Spatial model registry class python I/F   *
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
 * @file GModelSpatialRegistry.i
 * @brief GModelSpatialRegistry class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelSpatialRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelSpatialRegistry
 *
 * @brief Interface definition for the spatial model registry class.
 ***************************************************************************/
class GModelSpatialRegistry {
public:
    // Constructors and destructors
    GModelSpatialRegistry(void);
    GModelSpatialRegistry(const GModelSpatial* model);
    GModelSpatialRegistry(const GModelSpatialRegistry& registry);
    virtual ~GModelSpatialRegistry(void);

    // Methods
    int            size(void) const { return m_number; }
    GModelSpatial* alloc(const std::string& type) const;
    std::string    name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelSpatialRegistry class extension
 ***************************************************************************/
%extend GModelSpatialRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
