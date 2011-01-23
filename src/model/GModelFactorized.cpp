/***************************************************************************
 *            GModelFactorized.cpp  -  Model factorization class           *
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
 * @file GModelFactorized.cpp
 * @brief GModelFactorized class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GModelFactorized.hpp"
#include "GModelSpatialRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelTemporalRegistry.hpp"
#include "GModelTemporalConst.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_XML_SPATIAL           "GModelFactorized::xml_spatial(GXmlElement&)"
#define G_XML_SPECTRAL         "GModelFactorized::xml_spectral(GXmlElement&)"
#define G_XML_TEMPORAL         "GModelFactorized::xml_temporal(GXmlElement&)"

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
GModelFactorized::GModelFactorized(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model.
 ***************************************************************************/
GModelFactorized::GModelFactorized(const GModelFactorized& model)
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
GModelFactorized::~GModelFactorized(void)
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
 * @param[in] model Model.
 ***************************************************************************/
GModelFactorized& GModelFactorized::operator= (const GModelFactorized& model)
{
    // Execute only if object is not identical
    if (this != &model) {

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

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelFactorized::init_members(void)
{
    // Initialise members
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GModelFactorized::copy_members(const GModelFactorized& model)
{
    // Clone models
    m_spatial  = (model.m_spatial  != NULL) ? model.m_spatial->clone()  : NULL;
    m_spectral = (model.m_spectral != NULL) ? model.m_spectral->clone() : NULL;
    m_temporal = (model.m_temporal != NULL) ? model.m_temporal->clone() : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelFactorized::free_members(void)
{
    // Free memory
    if (m_spatial  != NULL) delete m_spatial;
    if (m_spectral != NULL) delete m_spectral;
    if (m_temporal != NULL) delete m_temporal;

    // Signal free pointers
    m_spatial  = NULL;
    m_spectral = NULL;
    m_temporal = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns list parameter pointers
 *
 * Gathers all parameter pointers from the model and returns them as a list
 * of GModelPar pointers.
 ***************************************************************************/
std::vector<GModelPar*> GModelFactorized::set_par_pointers(void)
{
    // Initialise parameters
    std::vector<GModelPar*> pars;

    // Determine number of parameters per type
    int n_spatial  = (spatial()  != NULL) ? spatial()->size()  : 0;
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spatial + n_spectral + n_temporal;

    // Continue only if there are parameters
    if (n_pars > 0) {

        // Gather spatial parameter pointers
        for (int i = 0; i < n_spatial; ++i)
            pars.push_back(&((*spatial())(i)));

        // Gather spectral parameters
        for (int i = 0; i < n_spectral; ++i)
            pars.push_back(&((*spectral())(i)));

        // Gather temporal parameters
        for (int i = 0; i < n_temporal; ++i)
            pars.push_back(&((*temporal())(i)));

    }

    // Return
    return pars;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the model information from the leafs "spectrum", "spatialModel",
 * and optionally "lightcurve". If the lightcurve leaf is not present a
 * constant model will be assumed. This is to assure compatibility with the
 * Fermi/LAT XML format.
 ***************************************************************************/
void GModelFactorized::xml_read(const GXmlElement& xml)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Initialise XML elements
    GXmlElement* spec = NULL;
    GXmlElement* spat = NULL;
    GXmlElement* temp = NULL;

    // Get pointers on spectrum and spatial model
    spec = (GXmlElement*)xml.element("spectrum", 0);
    spat = (GXmlElement*)xml.element("spatialModel", 0);

    // Clone spatial and spectral models
    m_spatial  = xml_spatial(*spat);
    m_spectral = xml_spectral(*spec);

    // Optionally get temporal model
    try {
        temp = (GXmlElement*)xml.element("lightcurve", 0);
        m_temporal = xml_temporal(*temp);
    }
    catch (GException::xml_name_not_found &e) {
        GModelTemporalConst temporal;
        m_temporal = temporal.clone();
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml Source library.
 * @param[in] name Model name.
 * @param[in] instruments Applicable instruments.
 *
 * Write the factorised model into XML document.
 ***************************************************************************/
void GModelFactorized::xml_write(GXmlElement& xml, const std::string& name,
                                 const std::string& instruments) const
{
    // Initialise pointer on source
    GXmlElement* src = NULL;

    // Search corresponding source
    int n = xml.elements("source");
    for (int k = 0; k < n; ++k) {
        GXmlElement* element = (GXmlElement*)xml.element("source", k);
        if (element->attribute("name") == name) {
            src = element;
            break;
        }
    }

    // If no source with corresponding name was found then append one
    if (src == NULL) {
        src = new GXmlElement("source");
        src->attribute("name") = name;
        if (spectral() != NULL) src->append(new GXmlElement("spectrum"));
        if (spatial()  != NULL) src->append(new GXmlElement("spatialModel"));
        xml.append(src);
    }

    // Set model type
    if (spatial() != NULL) {
        if (spatial()->type() == "PointSource")
            src->attribute("type", "PointSource");
        else
            src->attribute("type", "DiffuseSource");
    }
    else
        src->attribute("type", "unknown");

    // Set model attributes
    src->attribute("name", name);
    if (instruments.length() > 0)
        src->attribute("instrument", instruments);

    // Write spectral model
    if (spectral() != NULL) {
        GXmlElement* spec = (GXmlElement*)src->element("spectrum", 0);
        spectral()->write(*spec);
    }

    // Write spatial model
    if (spatial()) {
        GXmlElement* spat = (GXmlElement*)src->element("spatialModel", 0);
        spatial()->write(*spat);
    }

    // Write temporal model
    if (temporal()) {
        if (dynamic_cast<GModelTemporalConst*>(temporal()) == NULL) {
            GXmlElement* temp = (GXmlElement*)src->element("lightcurve", 0);
            temporal()->write(*temp);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Construct spatial model from XML element
 *
 * @param[in] spatial XML element containing spatial model information.
 *
 * @exception GException::model_invalid_spatial
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelSpatial* GModelFactorized::xml_spatial(const GXmlElement& spatial) const
{
    // Get spatial model type
    std::string type = spatial.attribute("type");

    // Get spatial model
    GModelSpatialRegistry registry;
    GModelSpatial*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL)
        ptr->read(spatial);

    // ... otherwise throw an exception
    else
        throw GException::model_invalid_spatial(G_XML_SPATIAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct spectral model from XML element
 *
 * @param[in] spectral XML element containing spectral model information.
 *
 * @exception GException::model_invalid_spectral
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelSpectral* GModelFactorized::xml_spectral(const GXmlElement& spectral) const
{
    // Get spectral model type
    std::string type = spectral.attribute("type");

    // Get spectral model
    GModelSpectralRegistry registry;
    GModelSpectral*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL)
        ptr->read(spectral);

    // ... otherwise throw an exception
    else
        throw GException::model_invalid_spectral(G_XML_SPECTRAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Construct temporal model from XML element
 *
 * @param[in] temporal XML element containing temporal model information.
 *
 * @exception GException::model_invalid_temporal
 *            Invalid spatial model type encountered.
 ***************************************************************************/
GModelTemporal* GModelFactorized::xml_temporal(const GXmlElement& temporal) const
{
    // Get temporal model type
    std::string type = temporal.attribute("type");

    // Get temporal model
    GModelTemporalRegistry registry;
    GModelTemporal*        ptr = registry.alloc(type);

    // If model if valid then read model from XML file
    if (ptr != NULL)
        ptr->read(temporal);

    // ... otherwise throw an exception
    else
        throw GException::model_invalid_temporal(G_XML_TEMPORAL, type);

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Verifies if model has all components
 ***************************************************************************/
bool GModelFactorized::valid_model(void) const
{
    // Set result
    bool result = ((m_spatial  != NULL) &&
                   (m_spectral != NULL) &&
                   (m_temporal != NULL));

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GModelFactorized::print_model(void) const
{
    // Initialise result string
    std::string result;

    // Determine number of parameters per type
    int n_spatial  = (spatial()  != NULL) ? spatial()->size()  : 0;
    int n_spectral = (spectral() != NULL) ? spectral()->size() : 0;
    int n_temporal = (temporal() != NULL) ? temporal()->size() : 0;
    int n_pars     = n_spatial + n_spectral + n_temporal;

    // Append model type
    result.append(parformat("Model type"));
    if (n_spatial > 0) {
        result.append("\""+spatial()->type()+"\"");
        if (n_spectral > 0 || n_temporal > 0)
            result.append(" * ");
    }
    if (n_spectral > 0) {
        result.append("\""+spectral()->type()+"\"");
        if (n_temporal > 0)
            result.append(" * ");
    }
    if (n_temporal > 0)
        result.append("\""+temporal()->type()+"\"");

    // Append parameters
    result.append("\n"+parformat("Number of parameters")+str(n_pars));
    result.append("\n"+parformat("Number of spatial par's")+str(n_spatial));
    for (int i = 0; i < n_spatial; ++i)
        result.append("\n"+(*spatial())(i).print());
    result.append("\n"+parformat("Number of spectral par's")+str(n_spectral));
    for (int i = 0; i < n_spectral; ++i)
        result.append("\n"+(*spectral())(i).print());
    result.append("\n"+parformat("Number of temporal par's")+str(n_temporal));
    for (int i = 0; i < n_temporal; ++i)
        result.append("\n"+(*temporal())(i).print());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                  Friends                                =
 =                                                                         =
 ==========================================================================*/
