/***************************************************************************
 *        GModelSpectralFunc.cpp  -  Spectral function model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralFunc.cpp
 * @brief GModelSpectralFunc class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralFunc.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                  "GModelSpectralFunc::par(int)"
#define G_READ                       "GModelSpectralFunc::read(GXmlElement&)"
#define G_WRITE                     "GModelSpectralFunc::write(GXmlElement&)"

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
GModelSpectralFunc::GModelSpectralFunc(void) : GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] filename File name of nodes.
 *
 * Construct an arbitrary spectral function from a list of nodes that is
 * found in the specified file.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const std::string& filename) :
                    GModelSpectral()
{
    // Initialise private members for clean destruction
    init_members();

    // Load nodes
    load_nodes(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element.
 *
 * Construct an arbitrary spectral function from a XML element.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const GXmlElement& xml) :
                    GModelSpectral()
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
 * @param[in] model Spectral component.
 ***************************************************************************/
GModelSpectralFunc::GModelSpectralFunc(const GModelSpectralFunc& model) :
                    GModelSpectral(model)
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
GModelSpectralFunc::~GModelSpectralFunc(void)
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
 * @param[in] model Spectral component.
 ***************************************************************************/
GModelSpectralFunc& GModelSpectralFunc::operator= (const GModelSpectralFunc& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

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
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GModelSpectralFunc::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpectralFunc* GModelSpectralFunc::clone(void) const
{
    return new GModelSpectralFunc(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 ***************************************************************************/
double GModelSpectralFunc::eval(const GEnergy& srcEng)
{
    // Interpolate function
    double func = m_nodes.interpolate(srcEng.MeV(), m_values);

    // Compute function value
    double value  = norm() * func;

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 *
 * The spectral model is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization of the function.
 * Note that the normalization is factorised into a scaling factor and a
 * value and that the method is expected to return the gradient with respect
 * to the parameter value \f$n_v\f$.
 *
 * The partial derivative of the normalization value is given by
 * \f[dI/dn_v=n_s f(E)\f]
 ***************************************************************************/
double GModelSpectralFunc::eval_gradients(const GEnergy& srcEng)
{
    // Interpolate function (dummy so far)
    double func = m_nodes.interpolate(srcEng.MeV(), m_values);

    // Compute function value
    double value  = norm() * func;

    // Compute partial derivatives of the parameter values
    double g_norm  = (m_norm.isfree())  ? m_norm.scale() * func : 0.0;

    // Set gradients
    m_norm.gradient(g_norm);

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing power law model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter name found in XML element.
 *
 * Read the function information from an XML element and load the nodes
 * from the associated file. The XML element is required to have an
 * attribute "file" that specifies the function nodes and a parameter
 * named "Normalization".
 ***************************************************************************/
void GModelSpectralFunc::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Spectral function requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Get value
    if (par->attribute("name") == "Normalization") {
        m_norm.read(*par);
    }
    else
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require \"Normalization\" parameter.");

    // Load nodes from file
    load_nodes(xml.attribute("file"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "FileFunction"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the spectral function information into an XML element. The XML
 * element has to be of type "FileFunction" and will have 1 parameter leaf
 * named "Normalization". Note that the function nodes will not be written
 * since they will not be altered by any method.
 ***************************************************************************/
void GModelSpectralFunc::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "FileFunction");

    // Verify model type
    if (xml.attribute("type") != "FileFunction")
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \"FileFunction\".");

    // If XML element has 0 nodes then append 1 parameter node
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"Normalization\""));
    }

    // Verify that XML element has exactly 1 parameter
    if (xml.elements() != 1 || xml.elements("parameter") != 1)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Spectral function requires exactly 1 parameter.");

    // Get parameter element
    GXmlElement* par = (GXmlElement*)xml.element("parameter", 0);

    // Set parameyter
    if (par->attribute("name") == "Normalization") {
        m_norm.write(*par);
    }
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Require \"Normalization\" parameter.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print powerlaw information
 ***************************************************************************/
std::string GModelSpectralFunc::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralFunc ===\n");
    result.append(parformat("Function file")+m_filename);
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
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralFunc::init_members(void)
{
    // Initialise parameters
    m_npars  = 1;
    m_par[0] = &m_norm;

    // Initialise powerlaw normalisation
    m_norm = GModelPar();
    m_norm.name("Normalization");
    m_norm.min(0.0);
    m_norm.max(1000.0);
    m_norm.value(1.0);
    m_norm.scale(1.0);
    m_norm.free();

    // Initialise other members
    m_nodes.clear();
    m_values.clear();
    m_filename.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral model component.
 ***************************************************************************/
void GModelSpectralFunc::copy_members(const GModelSpectralFunc& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_norm     = model.m_norm;
    m_nodes    = model.m_nodes;
    m_values   = model.m_values;
    m_filename = model.m_filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralFunc::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load nodes from file
 *
 * @param[in] filename File name.
 *
 * Node energies are stored in units of MeV.
 * At least 2 nodes are required.
 *
 * @todo Not yet implemented. 
 ***************************************************************************/
void GModelSpectralFunc::load_nodes(const std::string& filename)
{
    // Clear nodes and values
    m_nodes.clear();
    m_values.clear();

    // Set filename
    m_filename = filename;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] model Model.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpectralFunc& model)
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
GLog& operator<< (GLog& log, const GModelSpectralFunc& model)
{
    // Write spectrum into logger
    log << model.print();

    // Return logger
    return log;
}
