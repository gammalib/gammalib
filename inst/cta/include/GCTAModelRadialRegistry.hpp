/***************************************************************************
 *     GCTAModelRadialRegistry.hpp  -  CTA Radial model registry class     *
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
 * @file GCTAModelRadialRegistry.hpp
 * @brief GCTAModelRadialRegistry class interface definition
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRADIALREGISTRY_HPP
#define GCTAMODELRADIALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GCTAModelRadial.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialRegistry
 *
 * @brief Interface definition for the CTA radial model registry class.
 *
 * The registry class allows the registration of radial models for CTA that
 * are not necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of radial models that are available throughout all linked libraries. To
 * register a radial model it is sufficient to add
 *  const GCTAModelRadialXXX      g_cta_radial_XXX_seed;
 *  const GCTAModelRadialRegistry g_cta_radial_XXX_registry(&g_cta_radial_XXX_seed);
 * at the top of the .cpp file of the radial model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GCTAModelRadialRegistry {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GCTAModelRadialRegistry& registry);
    friend GLog&         operator<<(GLog& log, const GCTAModelRadialRegistry& registry);

public:
    // Constructors and destructors
    GCTAModelRadialRegistry(void);
    GCTAModelRadialRegistry(const GCTAModelRadial* model);
    GCTAModelRadialRegistry(const GCTAModelRadialRegistry& registry);
    virtual ~GCTAModelRadialRegistry(void);

    // Operators
    GCTAModelRadialRegistry& operator= (const GCTAModelRadialRegistry& registry);

    // Methods
    int              size(void) const { return m_number; }
    GCTAModelRadial* alloc(const std::string& type) const;
    std::string      name(const int& index) const;
    std::string      print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadialRegistry& registry);
    void free_members(void);

private:
    // Pricate members
    static int                     m_number;   //!< Number of models in registry
    static std::string*            m_names;    //!< Model names
    static const GCTAModelRadial** m_models;   //!< Pointer to seed models
};

#endif /* GCTAMODELRADIALREGISTRY_HPP */
