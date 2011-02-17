/***************************************************************************
 *     GModelRadialRegistry.hpp  -  Radial spatial model registry class    *
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
 * @file GModelRadialRegistry.hpp
 * @brief Radial spatial model class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELRADIALREGISTRY_HPP
#define GMODELRADIALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelRadial.hpp"


/***********************************************************************//**
 * @class GModelRadialRegistry
 *
 * @brief Radial spatial model registry class
 *
 * The registry class allows the registration of radial spatial models that
 * are not necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of radial spatial models that are available throughout all linked
 * libraries. To register a radial spatial model it is sufficient to add
 *  const GModelRadialXXX      g_radial_XXX_seed;
 *  const GModelRadialRegistry g_radial_XXX_registry(&g_radial_XXX_seed);
 * at the top of the .cpp file of the radial spatial model. Here, XXX is a
 * unique name that describes the model.
 ***************************************************************************/
class GModelRadialRegistry {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModelRadialRegistry& registry);
    friend GLog&         operator<<(GLog& log,        const GModelRadialRegistry& registry);

public:
    // Constructors and destructors
    GModelRadialRegistry(void);
    explicit GModelRadialRegistry(const GModelRadial* model);
    GModelRadialRegistry(const GModelRadialRegistry& registry);
    virtual ~GModelRadialRegistry(void);

    // Operators
    GModelRadialRegistry& operator=(const GModelRadialRegistry& registry);

    // Methods
    int           size(void) const { return m_number; }
    GModelRadial* alloc(const std::string& type) const;
    std::string   name(const int& index) const;
    std::string   print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelRadialRegistry& registry);
    void free_members(void);

private:
    // Pricate members
    static int                  m_number;   //!< Number of models in registry
    static std::string*         m_names;    //!< Model names
    static const GModelRadial** m_models;   //!< Pointer to seed models
};

#endif /* GMODELRADIALREGISTRY_HPP */
