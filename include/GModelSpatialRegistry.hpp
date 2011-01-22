/***************************************************************************
 *       GModelSpatialRegistry.hpp  -  Spatial model registry class        *
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
 * @file GModelSpatialRegistry.hpp
 * @brief GModelSpatialRegistry class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALREGISTRY_HPP
#define GMODELSPATIALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelSpatial.hpp"


/***********************************************************************//**
 * @class GModelSpatialRegistry
 *
 * @brief Interface definition for the spatial model registry class.
 *
 * The registry class allows the registration of spatial models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of spatial models that are available throughout all linked libraries. To
 * register a spatial model it is sufficient to add
 *  const GModelSpatialXXX      g_spatial_XXX_seed;
 *  const GModelSpatialRegistry g_spatial_XXX_registry(&g_spatial_XXX_seed);
 * at the top of the .cpp file of the spatial model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GModelSpatialRegistry {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModelSpatialRegistry& registry);
    friend GLog&         operator<<(GLog& log, const GModelSpatialRegistry& registry);

public:
    // Constructors and destructors
    GModelSpatialRegistry(void);
    GModelSpatialRegistry(const GModelSpatial* model);
    GModelSpatialRegistry(const GModelSpatialRegistry& registry);
    virtual ~GModelSpatialRegistry(void);

    // Operators
    GModelSpatialRegistry& operator= (const GModelSpatialRegistry& registry);

    // Methods
    int            size(void) const { return m_number; }
    GModelSpatial* alloc(const std::string& type) const;
    std::string    name(const int& index) const;
    std::string    print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialRegistry& registry);
    void free_members(void);

private:
    // Pricate members
    static int                   m_number;   //!< Number of models in registry
    static std::string*          m_names;    //!< Model names
    static const GModelSpatial** m_models;   //!< Pointer to seed models
};

#endif /* GMODELSPATIALREGISTRY_HPP */
