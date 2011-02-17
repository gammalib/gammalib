/***************************************************************************
 *            GModelPointSource.hpp  -  Point source model class           *
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
 * @file GModelPointSource.cpp
 * @brief GModelPointSource class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelRegistry.hpp"
#include "GModelPointSource.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelTemporalConst.hpp"

/* __ Globals ____________________________________________________________ */
const GModelPointSource g_pointsource_seed;
const GModelRegistry    g_pointsource_registry(&g_pointsource_seed);

/* __ Method name definitions ____________________________________________ */
#define G_DIR                                      "GModelPointSource::dir()"

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
GModelPointSource::GModelPointSource(void) : GModelSky()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GModelPointSource::GModelPointSource(const GXmlElement& xml) : GModelSky(xml)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct point source model from spectral component
 *
 * @param[in] ptsrc Point source model component.
 * @param[in] spectral Spectral model component.
 ***************************************************************************/
GModelPointSource::GModelPointSource(const GModelSpatialPtsrc& ptsrc, 
                                     const GModelSpectral& spectral) : 
                                     GModelSky()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_spatial  = ptsrc.clone();
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct point source model from spatial and spectral XML elements
 *
 * @param[in] ptsrc Point source XML element.
 * @param[in] spectral Spectral XML element.
 *
 * @todo Verify that spatial model is a point source.
 ***************************************************************************/
GModelPointSource::GModelPointSource(const GXmlElement& ptsrc, 
                                     const GXmlElement& spectral) : 
                                     GModelSky(ptsrc, spectral)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Point source model.
 ***************************************************************************/
GModelPointSource::GModelPointSource(const GModelPointSource& model) : 
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
GModelPointSource::~GModelPointSource(void)
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
 * @param[in] model Point source model.
 ***************************************************************************/
GModelPointSource& GModelPointSource::operator= (const GModelPointSource& model)
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
void GModelPointSource::clear(void)
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
GModelPointSource* GModelPointSource::clone(void) const
{
    return new GModelPointSource(*this);
}


/***********************************************************************//**
 * @brief Return source location
 *
 * @exception GException::no_point_source
 *            No point source model component found.
 ***************************************************************************/
GSkyDir GModelPointSource::dir(void) const
{
    // Get pointer on point source spatial model
    GModelSpatialPtsrc* ptr = dynamic_cast<GModelSpatialPtsrc*>(m_spatial);

    // Throw an exception if the spatial model is not a point source
    if (ptr == NULL)
        GException::no_point_source(G_DIR, name());

    // Return source location
    return (ptr->dir());
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GModelPointSource::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelPointSource ===");

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
void GModelPointSource::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Point source model.
 ***************************************************************************/
void GModelPointSource::copy_members(const GModelPointSource& model)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelPointSource::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
