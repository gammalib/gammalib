/***************************************************************************
 *          GModelDiffuseSource.cpp  -  Diffuse source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelDiffuseSource.cpp
 * @brief GModelDiffuseSource class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelRegistry.hpp"
#include "GModelDiffuseSource.hpp"
#include "GModelTemporalConst.hpp"

/* __ Globals ____________________________________________________________ */
const GModelDiffuseSource g_diffusesource_seed;
const GModelRegistry      g_diffusesource_registry(&g_diffusesource_seed);

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelDiffuseSource::GModelDiffuseSource(void) : GModelSky()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GModelDiffuseSource::GModelDiffuseSource(const GXmlElement& xml) : GModelSky(xml)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct diffuse source model from spatial and spectral components
 *
 * @param[in] spatial Spatial model.
 * @param[in] spectral Spectral model.
 ***************************************************************************/
GModelDiffuseSource::GModelDiffuseSource(const GModelSpatial& spatial,
                                         const GModelSpectral& spectral) : 
                                         GModelSky()
{
    // Initialise private members for clean destruction
    init_members();

    // Allocate temporal constant
    GModelTemporalConst temporal;

    // Clone model components
    m_spatial  = spatial.clone();
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct diffuse source model from spatial and spectral XML elements
 *
 * @param[in] spatial Spatial XML element.
 * @param[in] spectral Spectral XML element.
 ***************************************************************************/
GModelDiffuseSource::GModelDiffuseSource(const GXmlElement& spatial, 
                                         const GXmlElement& spectral) : 
                                         GModelSky(spatial, spectral)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Diffuse source model.
 ***************************************************************************/
GModelDiffuseSource::GModelDiffuseSource(const GModelDiffuseSource& model) : 
                                         GModelSky(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelDiffuseSource::~GModelDiffuseSource(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Diffuse source model.
 ***************************************************************************/
GModelDiffuseSource& GModelDiffuseSource::operator= (const GModelDiffuseSource& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSky::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the instance to an initial state.
 ***************************************************************************/
void GModelDiffuseSource::clear(void)
{
    // Free class members
    free_members();
    this->GModelSky::free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    this->GModelSky::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GModelDiffuseSource* GModelDiffuseSource::clone(void) const
{
    return new GModelDiffuseSource(*this);
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GModelDiffuseSource::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelDiffuseSource ===");

    // Append model
    result.append("\n"+print_model());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelDiffuseSource::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Diffuse source model.
 ***************************************************************************/
void GModelDiffuseSource::copy_members(const GModelDiffuseSource& model)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelDiffuseSource::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
