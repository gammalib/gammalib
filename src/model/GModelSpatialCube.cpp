/***************************************************************************
 *          GModelSpatialCube.cpp  -  Spatial map cube model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialCube.cpp
 * @brief GModelSpatialCube class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialCube.hpp"

/* __ Constants __________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_READ                        "GModelSpatialCube::read(GXmlElement&)"
#define G_WRITE                      "GModelSpatialCube::write(GXmlElement&)"

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
GModelSpatialCube::GModelSpatialCube(void) : GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 ***************************************************************************/
GModelSpatialCube::GModelSpatialCube(const GXmlElement& xml) :
                   GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spatial map cube component.
 ***************************************************************************/
GModelSpatialCube::GModelSpatialCube(const GModelSpatialCube& model) :
                   GModelSpatial(model)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpatialCube::~GModelSpatialCube(void)
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
 * @param[in] model Spatial map cube component.
 ***************************************************************************/
GModelSpatialCube& GModelSpatialCube::operator= (const GModelSpatialCube& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
void GModelSpatialCube::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpatialCube* GModelSpatialCube::clone(void) const
{
    return new GModelSpatialCube(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GModelSpatialCube::eval(const GSkyDir& srcDir)
{
    // Return value
    return 1.0;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GModelSpatialCube::eval_gradients(const GSkyDir& srcDir)
{
    // Set gradient to 0
    m_value.gradient(0.0);

    // Return value
    return 1.0;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the map cube information from an XML element. The XML element is
 * required to have 1 parameter named either "Normalization" or "Value".
 ***************************************************************************/
void GModelSpatialCube::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameters
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Map cube spatial model requires exactly 1 parameter.");

    // Get pointer on model parameter
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value")
        m_value.read(*par);
    else
        throw GException::model_invalid_parnames(G_READ, xml,
              "Map cube spatial model requires either \"Value\" or"
              " \"Normalization\" parameter.");

    // Load map cube from file
    load_cube(xml.attribute("file"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type "MapCubeFunction"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the map cube information into an XML element. The XML element has to
 * be of type "MapCubeFunction" and will have 1 parameter leaf named either
 * "Value" or "Normalization" (default).
 ***************************************************************************/
void GModelSpatialCube::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "MapCubeFunction");

    // Verify model type
    if (xml.attribute("type") != "MapCubeFunction")
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \"MapCubeFunction\".");

    // If XML element has 0 nodes then append parameter node. The name
    // of the node is "Normalization" as this is the Fermi-LAT standard.
    // We thus assure that the XML files will be compatible with
    // Fermi-LAT.
    if (xml.elements() == 0)
        xml.append(new GXmlElement("parameter name=\"Normalization\""));

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Map cube spatial model requires exactly 1 parameter.");

    // Get pointers on model parameter
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Set or update parameter
    if (par->attribute("name") == "Normalization" ||
        par->attribute("name") == "Value")
        m_value.write(*par);
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "Map cube spatial model requires either \"Value\" or"
              " \"Normalization\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print isotropic source information
 ***************************************************************************/
std::string GModelSpatialCube::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialCube ===\n");
    result.append(parformat("Map cube file")+m_filename);
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+parformat("Parameter "+str(i+1)));
        result.append(m_par[i]->print());
    }

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
void GModelSpatialCube::init_members(void)
{
    // Initialise parameters
    m_npars  = 1;
    m_par[0] = &m_value;

    // Initialise Value
    m_value = GModelPar();
    m_value.name("Normalization");
    m_value.value(1.0);
    m_value.scale(1.0);
    m_value.range(0.001, 1000.0);
    m_value.fix();

    // Initialise other members
    m_filename.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpatialCube members which should be copied.
 ***************************************************************************/
void GModelSpatialCube::copy_members(const GModelSpatialCube& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_value    = model.m_value;
    m_filename = model.m_filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialCube::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load map cube from file
 *
 * @param[in] filename File name.
 *
 * @todo Not yet implemented. 
 ***************************************************************************/
void GModelSpatialCube::load_cube(const std::string& filename)
{
    // Set filename
    m_filename = filename;

    // Return
    return;
}

/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] model Model.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpatialCube& model)
{
     // Write spectrum in output stream
    os << model.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] model Model.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GModelSpatialCube& model)
{
    // Write spectrum into logger
    log << model.print();

    // Return logger
    return log;
}
