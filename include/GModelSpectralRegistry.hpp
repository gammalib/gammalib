/***************************************************************************
 *      GModelSpectralRegistry.hpp  -  Spectral model registry class       *
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
 * @file GModelSpectralRegistry.hpp
 * @brief GModelSpectralRegistry class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRALREGISTRY_HPP
#define GMODELSPECTRALREGISTRY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GModelSpectral.hpp"


/***********************************************************************//**
 * @class GModelSpectralRegistry
 *
 * @brief Interface definition for the spectral model registry class.
 *
 * The registry class allows the registration of spectral models that are not
 * necessarily compiled into the GammaLib. It uses the static members
 * m_number, m_names, and m_models which are allocated globally to keep track
 * of spectral models that are available throughout all linked libraries. To
 * register a spectral model it is sufficient to add
 *  const GModelSpectralXXX      g_spectral_XXX_seed;
 *  const GModelSpectralRegistry g_spectral_XXX_registry(&g_spectral_XXX_seed);
 * at the top of the .cpp file of the spectral model. Here, XXX is a unique
 * name that describes the model.
 ***************************************************************************/
class GModelSpectralRegistry {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GModelSpectralRegistry& registry);
    friend GLog&         operator<<(GLog& log, const GModelSpectralRegistry& registry);

public:
    // Constructors and destructors
    GModelSpectralRegistry(void);
    GModelSpectralRegistry(const GModelSpectral* model);
    GModelSpectralRegistry(const GModelSpectralRegistry& registry);
    virtual ~GModelSpectralRegistry(void);

    // Operators
    GModelSpectralRegistry& operator= (const GModelSpectralRegistry& registry);

    // Methods
    int             size(void) const { return m_number; }
    GModelSpectral* alloc(const std::string& type) const;
    std::string     name(const int& index) const;
    std::string     print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralRegistry& registry);
    void free_members(void);

private:
    // Pricate members
    static int                    m_number;   //!< Number of models in registry
    static std::string*           m_names;    //!< Model names
    static const GModelSpectral** m_models;   //!< Pointer to seed models
};

#endif /* GMODELSPECTRALREGISTRY_HPP */
