/***************************************************************************
 *        GModelExtendedSource.cpp  -  Extended source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GModelExtendedSource.cpp
 * @brief Extended source model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelRegistry.hpp"
#include "GModelTemporalConst.hpp"

/* __ Globals ____________________________________________________________ */
const GModelExtendedSource g_extendedsource_seed;
const GModelRegistry       g_extendedsource_registry(&g_extendedsource_seed);

/* __ Method name definitions ____________________________________________ */
#define G_DIR                                   "GModelExtendedSource::dir()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelExtendedSource::GModelExtendedSource(void) : GModelSky()
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
GModelExtendedSource::GModelExtendedSource(const GXmlElement& xml) :
                      GModelSky(xml)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct extended source model from radial and spectral components
 *
 * @param[in] radial Radial model component.
 * @param[in] spectral Spectral model component.
 ***************************************************************************/
GModelExtendedSource::GModelExtendedSource(const GModelRadial&   radial,
                                           const GModelSpectral& spectral) :
                      GModelSky()
{
    // Initialise members
    init_members();

    // Allocate temporal constant model
    GModelTemporalConst temporal;

    // Clone model components
    m_spatial  = radial.clone();
    m_spectral = spectral.clone();
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct extended source model from radial and spectral XML
 *        elements
 *
 * @param[in] radial Radial XML element.
 * @param[in] spectral Spectral XML element.
 ***************************************************************************/
GModelExtendedSource::GModelExtendedSource(const GXmlElement& radial,
                                           const GXmlElement& spectral) :
                      GModelSky(radial, spectral)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Extended source model.
 ***************************************************************************/
GModelExtendedSource::GModelExtendedSource(const GModelExtendedSource& model) :
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
GModelExtendedSource::~GModelExtendedSource(void)
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
 * @param[in] model Extended source model.
 * @return Extended source model.
 ***************************************************************************/
GModelExtendedSource& GModelExtendedSource::operator=(const GModelExtendedSource& model)
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
void GModelExtendedSource::clear(void)
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
 *
 * @return Pointer to deep copy of extended source model.
 ***************************************************************************/
GModelExtendedSource* GModelExtendedSource::clone(void) const
{
    return new GModelExtendedSource(*this);
}


/***********************************************************************//**
 * @brief Print model information
 *
 * @return String containing model information.
 ***************************************************************************/
std::string GModelExtendedSource::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelExtendedSource ===");

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
void GModelExtendedSource::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Extended source model.
 ***************************************************************************/
void GModelExtendedSource::copy_members(const GModelExtendedSource& model)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelExtendedSource::free_members(void)
{
    // Return
    return;
}
