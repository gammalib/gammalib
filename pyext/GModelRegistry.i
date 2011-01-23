/***************************************************************************
 *           GModelRegistry.i  -  Model registry class python I/F          *
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
 * @file GModelRegistry.i
 * @brief GModelRegistry class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelRegistry.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelRegistry
 *
 * @brief Interface definition for the model registry class.
 ***************************************************************************/
class GModelRegistry {
public:
    // Constructors and destructors
    GModelRegistry(void);
    GModelRegistry(const GModel* model);
    GModelRegistry(const GModelRegistry& registry);
    virtual ~GModelRegistry(void);

    // Methods
    int         size(void) const { return m_number; }
    GModel*     alloc(const std::string& type) const;
    std::string name(const int& index) const;
};


/***********************************************************************//**
 * @brief GModelRegistry class extension
 ***************************************************************************/
%extend GModelRegistry {
    char *__str__() {
        return tochar(self->print());
    }
};
