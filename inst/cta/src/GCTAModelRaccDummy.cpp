/***************************************************************************
 *     GCTAModelRaccDummy.cpp  -  Radial acceptance dummy model class      *
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
 * @file GCTAModelRaccDummy.cpp
 * @brief GCTAModelRaccDummy class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelRegistry.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAModelRaccDummy.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelRaccDummy g_cta_racc_dummy_seed;
const GModelRegistry     g_cta_racc_dummy_registry(&g_cta_racc_dummy_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                       "GCTAModelRaccDummy::read(GXmlElement&)"
#define G_WRITE                     "GCTAModelRaccDummy::write(GXmlElement&)"

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
GCTAModelRaccDummy::GCTAModelRaccDummy(void) : GModelData(),
                                               GModelFactorized()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GCTAModelRaccDummy::GCTAModelRaccDummy(const GXmlElement& xml) :
                                       GModelData(xml),
                                       GModelFactorized()
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
GCTAModelRaccDummy::GCTAModelRaccDummy(const GCTAModelRaccDummy& model) :
                                       GModelData(model),
                                       GModelFactorized(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAModelRaccDummy::~GCTAModelRaccDummy(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Radial acceptance model.
 ***************************************************************************/
GCTAModelRaccDummy& GCTAModelRaccDummy::operator= (const GCTAModelRaccDummy& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelData::operator=(model);
        this->GModelFactorized::operator=(model);

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
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GCTAModelRaccDummy::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelFactorized::free_members();
    this->GModelData::free_members();
    this->GModel::free_members();

    // Initialise members
    this->GModel::init_members();
    this->GModelData::init_members();
    this->GModelFactorized::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAModelRaccDummy* GCTAModelRaccDummy::clone(void) const
{
    return new GCTAModelRaccDummy(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * @todo Verify that CTA instrument direction pointer is valid
 ***************************************************************************/
double GCTAModelRaccDummy::eval(const GEvent& event, const GObservation& obs)
{
    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Evaluate function
    double value = 1.0;
    if (spatial()  != NULL) value *= spatial()->eval(cta_dir->skydir());
    if (spectral() != NULL) value *= spectral()->eval(event.energy());
    if (temporal() != NULL) value *= temporal()->eval(event.time());

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] event Observed event.
 * @param[in] obs Observation.
 *
 * @todo Verify that CTA instrument direction pointer is valid
 ***************************************************************************/
double GCTAModelRaccDummy::eval_gradients(const GEvent& event,
                                          const GObservation& obs)
{
    // Get instrument direction
    const GInstDir*    inst_dir = &(event.dir());
    const GCTAInstDir* cta_dir  = dynamic_cast<const GCTAInstDir*>(inst_dir);

    // Evaluate function and gradients
    double value = 1.0;
    if (spatial()  != NULL) value *= spatial()->eval_gradients(cta_dir->skydir());
    if (spectral() != NULL) value *= spectral()->eval_gradients(event.energy());
    if (temporal() != NULL) value *= temporal()->eval_gradients(event.time());

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GCTAModelRaccDummy::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Read model
    xml_read(xml);

    // Store name and instruments
    this->name(xml.attribute("name"));
    this->instruments(xml.attribute("instrument"));

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 ***************************************************************************/
void GCTAModelRaccDummy::write(GXmlElement& xml) const
{
    // Write model
    xml_write(xml, name(), instruments());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print model information
 ***************************************************************************/
std::string GCTAModelRaccDummy::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAModelRaccDummy ===");

    // Append name and instrument
    result.append("\n"+print_name_instrument());

    // Append model
    result.append("\n"+print_model());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelRaccDummy::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Model.
 ***************************************************************************/
void GCTAModelRaccDummy::copy_members(const GCTAModelRaccDummy& model)
{
    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelRaccDummy::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * @todo Not yet implemented
 ***************************************************************************/
void GCTAModelRaccDummy::set_pointers(void)
{
    // Set pointers
    m_pars = set_par_pointers();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
