/***************************************************************************
 *        GModelExtendedSource.hpp  -  Extended source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelRegistry.hpp"
#include "GModelRadialRegistry.hpp"
#include "GModelTemporalConst.hpp"

/* __ Globals ____________________________________________________________ */
const GModelExtendedSource g_extendedsource_seed;
const GModelRegistry       g_extendedsource_registry(&g_extendedsource_seed);

/* __ Method name definitions ____________________________________________ */
#define G_DIR                                   "GModelExtendedSource::dir()"
#define G_XML_RADIAL         "GModelExtendedSource::xml_radial(GXmlElement&)"

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
GModelExtendedSource::GModelExtendedSource(const GXmlElement& xml)
                                           : GModelSky()
{
    // Initialise members
    init_members();

    // Get pointers on spectrum and spatial model
    GXmlElement* spec = static_cast<GXmlElement*>(xml.element("spectrum", 0));
    GXmlElement* rad  = static_cast<GXmlElement*>(xml.element("spatialModel", 0));

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_radial(*rad);
    m_spectral = xml_spectral(*spec);
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

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
                                           const GModelSpectral& spectral)
                                           : GModelSky()
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
                                           const GXmlElement& spectral)
                                           : GModelSky(radial, spectral)
{
    // Initialise members
    init_members();

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_radial(radial);
    m_spectral = xml_spectral(spectral);
    m_temporal = temporal.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Extended source model.
 ***************************************************************************/
GModelExtendedSource::GModelExtendedSource(const GModelExtendedSource& model)
                                           : GModelSky(model)
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
 ***************************************************************************/
GModelExtendedSource* GModelExtendedSource::clone(void) const
{
    return new GModelExtendedSource(*this);
}


/***********************************************************************//**
 * @brief Read extended model from XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GModelExtendedSource::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Get pointers on spectrum and spatial model
    GXmlElement* spec = static_cast<GXmlElement*>(xml.element("spectrum", 0));
    GXmlElement* rad  = static_cast<GXmlElement*>(xml.element("spatialModel", 0));

    // Allocate constant
    GModelTemporalConst temporal;

    // Clone spatial and spectral models
    m_spatial  = xml_radial(*rad);
    m_spectral = xml_spectral(*spec);
    m_temporal = temporal.clone();

    // Set model name
    name(xml.attribute("name"));

    // Set instruments
    instruments(xml.attribute("instrument"));

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write extended model into XML element
 *
 * @param[in] xml Source library.
 ***************************************************************************/
void GModelExtendedSource::write(GXmlElement& xml) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = static_cast<GXmlElement*>(xml.element("source", k));
        if (element->attribute("name") == name()) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = new GXmlElement("source");
        src->attribute("name") = name();
        if (spectral() != NULL) src->append(new GXmlElement("spectrum"));
        if (radial()   != NULL) src->append(new GXmlElement("spatialModel"));
        xml.append(src);
    }

    // Set model attributes
    src->attribute("name", name());
    src->attribute("type", type());
    std::string instruments = this->instruments();
    if (instruments.length() > 0) {
        src->attribute("instrument", instruments);
    }

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = static_cast<GXmlElement*>(src->element("spectrum", 0));
        spectral()->write(*spec);
    }

    // Write spatial model
    if (radial() != NULL) {
        GXmlElement* rad = static_cast<GXmlElement*>(src->element("spatialModel", 0));
        radial()->write(*rad);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer to radial component
 *
 * This method returns NULL is the spatial component is not a radial model.
 ***************************************************************************/
GModelRadial* GModelExtendedSource::radial(void) const
{
    // Get pointer on radial component
    GModelRadial* ptr = dynamic_cast<GModelRadial*>(m_spatial);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return source location
 *
 * @exception GException::no_extended_source
 *            No extended source model component found.
 ***************************************************************************/
GSkyDir GModelExtendedSource::dir(void) const
{
    // Get pointer on point radial model
    GModelRadial* ptr = dynamic_cast<GModelRadial*>(m_spatial);

    // Throw an exception if the spatial model is not a point source
    if (ptr == NULL) {
        GException::no_extended_source(G_DIR, name());
    }

    // Return source location
    return (ptr->dir());
}


/***********************************************************************//**
 * @brief Print model information
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
 * @param[in] model Point source model.
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


/***********************************************************************//**
 * @brief Construct radial model from XML element
 *
 * @param[in] radial XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Invalid radial spatial model type encountered.
 ***************************************************************************/
GModelRadial* GModelExtendedSource::xml_radial(const GXmlElement& radial) const
{
    // Get radial spatial model type
    std::string type = radial.attribute("type");

    // Get radial spatial model
    GModelRadialRegistry registry;
    GModelRadial*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL) {
        ptr->read(radial);
    }

    // ... otherwise throw an exception
    else {
        throw GException::model_invalid_spatial(G_XML_RADIAL, type);
    }

    // Return pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
